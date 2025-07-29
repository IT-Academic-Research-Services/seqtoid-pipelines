"""
reference-taxa-download.py: Downloads reference genome files from NCBI based on the names.dmp file
format as found in NCBI taxonomy.
That file can be found by the following.
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
Requires the NCBI datasets program from: https://www.ncbi.nlm.nih.gov/datasets/
Prefers genomes with the--refrence argument, falls back to --assembly-level chromosome,complete
"""

import os
import subprocess
import logging
import multiprocessing as mp
from pathlib import Path
from typing import Dict, Tuple, List, Optional
from tqdm import tqdm

SCIENTIFIC_NAME_TAG = 'scientific name'
SYNONYM_NAME_TAG = 'synonym'
LOG_FILE = os.path.join(os.getcwd(), "process_taxonomy_download_genomes.log")

logging.basicConfig(
    filename=LOG_FILE,
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)

def load_taxids_by_block(names_file):
    taxid_blocks = {}
    taxids_found = set()

    if not os.path.isfile(names_file):
        print(f"The names file {names_file} in .dmp format from NCBI taxaonomy cannot be found ")
        return None
    with open(names_file) as f:
        for fline in f:
            ff = fline.strip().split("\t|\t")
            if len(ff) >= 4:
                taxid = ff[0]

                taxids_found.add(taxid)
                name = ff[1].strip()
                name_type_base = ff[3].strip()
                name_type_no_pipe = name_type_base.strip().rstrip("|")
                name_type_make_sure = name_type_no_pipe.strip() # ugh

                if taxid in taxid_blocks:
                    taxid_blocks[taxid].append([name, name_type_make_sure])
                else:
                    taxid_blocks[taxid] = [[name, name_type_make_sure]]


    print(f"Total taxids in file {len(taxids_found)}")

    return taxid_blocks


def filter_taxid_blocks(taxid_blocks_dict):
    taxids = {}
    for taxid, name_data in taxid_blocks_dict.items():
        final_name = "N/A"
        final_name_type = "N/A"

        for name_datum in name_data:
            name_candidate = name_datum[0]
            name_type_candidate = name_datum[1]
            if SCIENTIFIC_NAME_TAG in name_type_candidate:

                final_name = name_candidate
                final_name_type = name_type_candidate
                break
            elif SYNONYM_NAME_TAG in name_type_candidate:
                final_name = name_candidate
                final_name_type = name_type_candidate
        taxids[taxid] = [final_name, final_name_type]
    return taxids




def download_genome_counter(args):
    taxid, name_data = args
    return download_genome(taxid, name_data)
def download_genome(taxid, name_data):
    cwd = os.getcwd()

    print(f"{taxid}")
    taxon_name = name_data[0]
    try:

        zip_file = f"ncbi_dataset_{taxid}.zip"
        cmd = ["datasets", "download", "genome", "taxon", taxid, "--reference", "--filename", zip_file]
        result = subprocess.run(
            cmd, cwd=cwd, capture_output=True, text=True, check=False
        )

        if result.returncode == 0:
            unzip_cmd = ["unzip", "-o", zip_file, "-d", f"genome_{taxid}"]
            unzip_result = subprocess.run(
                unzip_cmd, cwd=cwd, capture_output=True, text=True, check=False
            )
            if unzip_result.returncode == 0:
                genome_dir = os.path.join(cwd, f"genome_{taxid}", "ncbi_dataset", "data")
                for root, _, files in os.walk(genome_dir):
                    for file in files:
                        if file.endswith("_genomic.fna"):
                            genome_file = os.path.join(root, file)
                            new_file = os.path.join(cwd, f"taxid_{taxid}_{taxon_name.replace(' ', '_')}_genomic.fna")
                            os.rename(genome_file, new_file)
                            subprocess.run(["rm", "-rf", f"genome_{taxid}", zip_file], cwd=cwd)
                            logging.info(f"Successfully downloaded reference genome for TaxID {taxid}: {new_file}")
                        return taxid, taxon_name, new_file
                logging.warning(f"No genomic.fna file found for TaxID {taxid} in reference download")
                subprocess.run(["rm", "-rf", f"genome_{taxid}", zip_file], cwd=cwd)

        logging.info(f"Attempting fallback download for TaxID {taxid} (non-reference)")
        cmd_fallback = ["datasets", "download", "genome", "taxon", taxid, "--assembly-level", "chromosome,complete", "--filename", zip_file]
        result_fallback = subprocess.run(
            cmd_fallback, cwd=cwd, capture_output=True, text=True, check=False
        )

        if result_fallback.returncode == 0:
            genome_dir = os.path.join(cwd, f"genome_{taxid}", "ncbi_dataset", "data")
        for root, _, files in os.walk(cwd):
            for file in files:
                if file.endswith("_genomic.fna"):
                    genome_file = os.path.join(root, file)
                    new_file = os.path.join(cwd, f"taxid_{taxid}_{taxon_name.replace(' ', '_')}_genomic.fna")
                    os.rename(genome_file, new_file)
                    subprocess.run(["rm", "-rf", f"genome_{taxid}", zip_file], cwd=cwd)
                    logging.info(f"Successfully downloaded reference genome for TaxID {taxid}: {new_file}")
                return taxid, taxon_name, new_file
        logging.warning(f"No genomic.fna file found for TaxID {taxid} in reference download")
        subprocess.run(["rm", "-rf", f"genome_{taxid}", zip_file], cwd=cwd)

        logging.warning(f"No genome available for TaxID {taxid}")
        return taxid, taxon_name, "N/A"
    except Exception as e:
        logging.error(f"Error processing TaxID {taxid}: {e}")
        return taxid, taxon_name, "N/A"

def write_tsv(results: List[Tuple[str, str, str]], output_tsv) -> None:
    try:
        with open(output_tsv, 'w', encoding='utf-8') as f:
            f.write("TaxID\tScientificName\tGenomeFile\n")
            for taxid, name, genome_file in results:
                f.write(f"{taxid}\t{name}\t{genome_file}\n")
        logging.info(f"Wrote results to {output_tsv}")
    except Exception as e:
        logging.error(f"Error writing TSV file {output_tsv}: {e}")
        raise

def main(namesfile = 'names.dmp'):
    print("NCBI Reference Genome download")

    taxid_blocks = load_taxids_by_block(names_file=namesfile)


    taxids = filter_taxid_blocks(taxid_blocks)

    with open("filtered_taxids.tsv", 'w') as f:
        for taxid, name_data in taxids.items():
            f.write(str(taxid))
            f.write("\t")
            f.write(str(name_data[0]))
            f.write("\t")
            f.write(str(name_data[1]))
            f.write("\n")

    num_cores = mp.cpu_count()
    print(f"Using {num_cores} threads to process and download {len(taxids)} taxa.")


    with mp.Pool(processes=num_cores) as pool:
        results = list(tqdm(pool.imap_unordered(download_genome_counter, taxids.items()), total=len(taxids), desc="Processing taxon IDs"))

    # for taxid, name_data in taxids.items():
    #     taxid, taxon_name, new_file = download_genome(taxid, name_data)
    #     print(f"{taxid}    {taxon_name}    {new_file}")



if __name__ == "__main__":
    # main("names.dmp")
    main("names-head.dmp")
