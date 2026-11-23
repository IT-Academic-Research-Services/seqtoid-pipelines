#include <cctype>
#include <climits>
#include <cuda_runtime.h>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
#define CK(x) do{cudaError_t e=(x);if(e!=cudaSuccess)throw std::runtime_error(cudaGetErrorString(e));}while(0)
using Clock=std::chrono::steady_clock;
double ms(Clock::time_point t){return std::chrono::duration<double,std::milli>(Clock::now()-t).count();}
__constant__ int mat[576];
const int hostmat[576]={4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0,-2,-1,0,-4,-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3,-1,0,-1,-4,-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3,3,0,-1,-4,-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3,4,1,-1,-4,0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2,-4,-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2,0,3,-1,-4,-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2,1,4,-1,-4,0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3,-1,-2,-1,-4,-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3,0,0,-1,-4,-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3,-3,-3,-1,-4,-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1,-4,-3,-1,-4,-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2,0,1,-1,-4,-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1,-3,-1,-1,-4,-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1,-3,-3,-1,-4,-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2,-2,-1,-2,-4,1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2,0,0,0,-4,0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0,-1,-1,0,-4,-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3,-4,-3,-2,-4,-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1,-3,-2,-1,-4,0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4,-3,-2,-1,-4,-2,-1,3,4,-3,0,1,-1,0,-3,-4,0,-3,-3,-2,0,-1,-4,-3,-3,4,1,-1,-4,-1,0,0,1,-3,3,4,-2,0,-3,-3,1,-1,-3,-1,0,-1,-3,-2,-2,1,4,-1,-4,0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2,0,0,-2,-1,-1,-1,-1,-1,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,1};
const std::string alphabet="ARNDCQEGHILKMFPSTWYVBZX*";
struct Tile{size_t first,n,begin,end;};
__host__ __device__ unsigned long long mix(unsigned long long x){x+=0x9e3779b97f4a7c15ULL;x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL;x=(x^(x>>27))*0x94d049bb133111ebULL;return x^(x>>31);}
// One warp per protein; lanes scan independent ungapped diagonals.
__global__ void score(const unsigned char* ref,const uint64_t* off,const int* len,int n,const unsigned char* query,int qlen,int* scores){
 int lane=threadIdx.x%32;int t=(blockIdx.x*blockDim.x+threadIdx.x)/32;if(t>=n)return;
 int best=0,L=len[t];const unsigned char* s=ref+off[t];
 for(long long d=-(long long)qlen+1+lane;d<L;d+=32){int qi=d<0?(int)-d:0;int ti=d>0?(int)d:0;int run=0;
  for(;qi<qlen&&ti<L;++qi,++ti){run=max(0,run+mat[query[qi]*24+s[ti]]);best=max(best,run);}}
 for(int delta=16;delta;delta/=2)best=max(best,__shfl_down_sync(0xffffffff,best,delta));if(lane==0)scores[t]=best;
}
// One block per query/tile reduction: checksum covers every target score.
__global__ void reduce(const int* scores,int n,uint64_t first,unsigned long long* sums,int* bests,unsigned long long* ids,int q){
 __shared__ unsigned long long shsum[256],shid[256];__shared__ int shbest[256];int k=threadIdx.x;unsigned long long sum=0,id=~0ULL;int best=-1;
 for(int i=k;i<n;i+=256){int s=scores[i];unsigned long long key=first+i;sum+=mix(key)^mix((unsigned long long)s);if(s>best||(s==best&&key<id)){best=s;id=key;}}
 shsum[k]=sum;shbest[k]=best;shid[k]=id;__syncthreads();
 for(int d=128;d;d/=2){if(k<d){shsum[k]+=shsum[k+d];if(shbest[k+d]>shbest[k]||(shbest[k+d]==shbest[k]&&shid[k+d]<shid[k])){shbest[k]=shbest[k+d];shid[k]=shid[k+d];}}__syncthreads();}
 if(k==0){sums[q]+=shsum[0];if(shbest[0]>bests[q]||(shbest[0]==bests[q]&&shid[0]<ids[q])){bests[q]=shbest[0];ids[q]=shid[0];}}
}
int cpu(const unsigned char* s,int L,const unsigned char* q,int Q){int best=0;for(long long d=-(long long)Q+1;d<L;++d){int qi=d<0?-d:0,ti=d>0?d:0,run=0;for(;qi<Q&&ti<L;++qi,++ti){run=std::max(0,run+hostmat[q[qi]*24+s[ti]]);best=std::max(best,run);}}return best;}
int main(int argc,char**argv){try{
 std::string fasta,schedule="tile";int Q=32,qlen=31,B=8,device=0,seed=17,validate=0;size_t tilebytes=64ULL<<20,maxrecords=0;
 for(int i=1;i<argc;i++){std::string a=argv[i];if(i+1>=argc)throw std::runtime_error("Missing option value");std::string v=argv[++i];if(a=="--fasta")fasta=v;else if(a=="--schedule")schedule=v;else if(a=="--queries")Q=std::stoi(v);else if(a=="--qlen")qlen=std::stoi(v);else if(a=="--batch")B=std::stoi(v);else if(a=="--tile-mib")tilebytes=std::stoull(v)<<20;else if(a=="--max-records")maxrecords=std::stoull(v);else if(a=="--device")device=std::stoi(v);else if(a=="--seed")seed=std::stoi(v);else if(a=="--validate")validate=std::stoi(v);else throw std::runtime_error("Unknown option "+a);}
 if(Q<1||qlen<1||qlen>100000||B<1||tilebytes==0||(schedule!="query"&&schedule!="tile"&&schedule!="resident"))throw std::runtime_error("Invalid configuration");
 auto prep=Clock::now();std::ifstream in(fasta);if(!in)throw std::runtime_error("Cannot open FASTA");std::vector<unsigned char> ref;std::vector<uint64_t> offsets;std::vector<int> lengths;std::string line,seq;bool header=false;
 auto append=[&](){if(seq.empty())throw std::runtime_error("Empty protein");if(seq.size()>100000)throw std::runtime_error("Protein exceeds 100000 aa");offsets.push_back(ref.size());lengths.push_back(seq.size());for(unsigned char c:seq){size_t p=alphabet.find(std::toupper(c));ref.push_back(p==std::string::npos?22:p);}seq.clear();};
 while(std::getline(in,line)){if(!line.empty()&&line[0]=='>'){if(header){append();if(maxrecords&&lengths.size()>=maxrecords){header=false;break;}}header=true;}else{for(unsigned char c:line)if(!std::isspace(c)){if(!header)throw std::runtime_error("Invalid FASTA");seq+=c;}}}if(header)append();if(lengths.empty())throw std::runtime_error("No proteins");
 std::vector<size_t> eligible;for(size_t i=0;i<lengths.size();i++)if(lengths[i]>=qlen)eligible.push_back(i);if(eligible.empty())throw std::runtime_error("No target long enough for queries");std::mt19937 rng(seed);std::vector<unsigned char> queries((size_t)Q*qlen);for(int q=0;q<Q;q++){size_t t=eligible[rng()%eligible.size()];size_t start=rng()%(lengths[t]-qlen+1);std::copy_n(ref.data()+offsets[t]+start,qlen,queries.data()+(size_t)q*qlen);}
 std::vector<Tile> tiles;for(size_t i=0;i<lengths.size();){size_t j=i,begin=offsets[i],end=begin;while(j<lengths.size()){size_t e=offsets[j]+lengths[j];if(j>i&&e-begin+(j-i+1)*12>tilebytes)break;end=e;++j;}tiles.push_back({i,j-i,begin,end});i=j;}
 if(schedule=="resident")tiles={{0,lengths.size(),0,ref.size()}};
 size_t cap=0,ncap=0;for(auto t:tiles){cap=std::max(cap,t.end-t.begin);ncap=std::max(ncap,t.n);}if(ncap>INT32_MAX)throw std::runtime_error("Tile too many targets");double prepms=ms(prep);
 auto init=Clock::now();CK(cudaSetDevice(device));CK(cudaMemcpyToSymbol(mat,hostmat,sizeof(hostmat)));cudaDeviceProp prop;CK(cudaGetDeviceProperties(&prop,device));
 unsigned char *dr,*dq,*hr;uint64_t *dof,*hof;int *dl,*hl,*ds,*db;unsigned long long *dc,*di;
 CK(cudaMallocHost(&hr,cap));CK(cudaMallocHost(&hof,ncap*8));CK(cudaMallocHost(&hl,ncap*4));CK(cudaMalloc(&dr,cap));CK(cudaMalloc(&dof,ncap*8));CK(cudaMalloc(&dl,ncap*4));CK(cudaMalloc(&ds,ncap*4));CK(cudaMalloc(&dq,queries.size()));CK(cudaMalloc(&dc,Q*8));CK(cudaMalloc(&db,Q*4));CK(cudaMalloc(&di,Q*8));
 cudaEvent_t e0,e1,e2,e3;CK(cudaEventCreate(&e0));CK(cudaEventCreate(&e1));CK(cudaEventCreate(&e2));CK(cudaEventCreate(&e3));CK(cudaMemcpy(dq,queries.data(),queries.size(),cudaMemcpyHostToDevice));double initms=ms(init);
 auto upload=[&](Tile t){std::memcpy(hr,ref.data()+t.begin,t.end-t.begin);for(size_t i=0;i<t.n;i++){hof[i]=offsets[t.first+i]-t.begin;hl[i]=lengths[t.first+i];}};
 auto copy=[&](Tile t){CK(cudaMemcpyAsync(dr,hr,t.end-t.begin,cudaMemcpyHostToDevice));CK(cudaMemcpyAsync(dof,hof,t.n*8,cudaMemcpyHostToDevice));CK(cudaMemcpyAsync(dl,hl,t.n*4,cudaMemcpyHostToDevice));};
 auto launch=[&](Tile t,int q){score<<<(t.n+7)/8,256>>>(dr,dof,dl,(int)t.n,dq+(size_t)q*qlen,qlen,ds);CK(cudaGetLastError());};
 // Warm up and validate exact per-target scores, outside measured execution.
 Tile first=tiles.front();upload(first);copy(first);launch(first,0);CK(cudaDeviceSynchronize());
 if(validate){for(auto t:tiles){upload(t);copy(t);std::vector<int> got(t.n);for(int q=0;q<Q;q++){launch(t,q);CK(cudaMemcpy(got.data(),ds,t.n*4,cudaMemcpyDeviceToHost));for(size_t i=0;i<t.n;i++)if(got[i]!=cpu(ref.data()+offsets[t.first+i],lengths[t.first+i],queries.data()+(size_t)q*qlen,qlen))throw std::runtime_error("CPU/GPU score mismatch");}}}
 CK(cudaMemset(dc,0,Q*8));CK(cudaMemset(db,0xff,Q*4));CK(cudaMemset(di,0xff,Q*8));CK(cudaDeviceSynchronize());
 double staging=0,h2d=0,kernel=0,reduction=0;uint64_t bytes=0,uploads=0;auto wall=Clock::now();
 auto load=[&](Tile t){auto st=Clock::now();upload(t);staging+=ms(st);CK(cudaEventRecord(e0));copy(t);CK(cudaEventRecord(e1));CK(cudaEventSynchronize(e1));float v;CK(cudaEventElapsedTime(&v,e0,e1));h2d+=v;bytes+=t.end-t.begin+t.n*12;uploads++;};
 auto work=[&](Tile t,int q){CK(cudaEventRecord(e1));launch(t,q);CK(cudaEventRecord(e2));reduce<<<1,256>>>(ds,(int)t.n,t.first,dc,db,di,q);CK(cudaGetLastError());CK(cudaEventRecord(e3));CK(cudaEventSynchronize(e3));float v;CK(cudaEventElapsedTime(&v,e1,e2));kernel+=v;CK(cudaEventElapsedTime(&v,e2,e3));reduction+=v;};
 if(schedule=="query"){for(int q=0;q<Q;q++)for(auto t:tiles){load(t);work(t,q);}}
 else if(schedule=="tile"){for(int q=0;q<Q;q+=B)for(auto t:tiles){load(t);for(int j=q;j<std::min(Q,q+B);j++)work(t,j);}}
 else{load(tiles[0]);for(int q=0;q<Q;q++)work(tiles[0],q);}
 std::vector<unsigned long long> sums(Q),ids(Q);std::vector<int>bests(Q);auto dt=Clock::now();CK(cudaMemcpy(sums.data(),dc,Q*8,cudaMemcpyDeviceToHost));CK(cudaMemcpy(ids.data(),di,Q*8,cudaMemcpyDeviceToHost));CK(cudaMemcpy(bests.data(),db,Q*4,cudaMemcpyDeviceToHost));CK(cudaDeviceSynchronize());double d2h=ms(dt),wallms=ms(wall);
 std::cout<<"{\"schedule\":\""<<schedule<<"\",\"gpu\":\""<<prop.name<<"\",\"queries\":"<<Q<<",\"qlen\":"<<qlen<<",\"batch\":"<<B<<",\"targets\":"<<lengths.size()<<",\"residues\":"<<ref.size()<<",\"tiles\":"<<tiles.size()<<",\"tile_mib\":"<<(tilebytes>>20)<<",\"validated\":"<<(validate?"true":"false")<<",\"prep_ms\":"<<prepms<<",\"init_query_upload_ms\":"<<initms<<",\"staging_ms\":"<<staging<<",\"h2d_ms\":"<<h2d<<",\"kernel_ms\":"<<kernel<<",\"reduction_ms\":"<<reduction<<",\"d2h_wall_ms\":"<<d2h<<",\"wall_ms\":"<<wallms<<",\"reference_h2d_bytes\":"<<bytes<<",\"uploads\":"<<uploads<<",\"results\":[";
 for(int q=0;q<Q;q++){if(q)std::cout<<",";std::cout<<"["<<bests[q]<<","<<ids[q]<<",\""<<sums[q]<<"\"]";}std::cout<<"]}\n";
 return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}}
