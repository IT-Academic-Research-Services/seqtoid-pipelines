NVCC ?= nvcc
NVCCFLAGS ?= -O3 -std=c++17 -lineinfo -arch=sm_89

bench: bench.cu
	$(NVCC) $(NVCCFLAGS) $< -o $@

clean:
	rm -f bench
