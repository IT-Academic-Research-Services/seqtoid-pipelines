### Create DB (`--module create_db`)
This workflow creates or updates an HDF5 database for reference sequences, enabling fast retrieval by accession or sequence ID.

**Steps**:
- Parses input FASTA/FASTQ and indexes sequences.
- Stores in HDF5 format for efficient querying.

**Inputs**: Reference files or directories.
**Outputs**: HDF5 file with indexed sequences.