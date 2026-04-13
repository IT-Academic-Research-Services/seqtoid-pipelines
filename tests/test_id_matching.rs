use seqtoid_pipelines::utils::fastx::compare_read_ids;

#[test]
fn test_sra_id_matching() {
    assert!(compare_read_ids("SRR8073913.1.1", "SRR8073913.1.2"));
    
    // SRA-style headers with space and read number (often seen in fastq-dump)
    // @SRR8073913.1 1 length=101
    // @SRR8073913.1 2 length=101
    assert!(compare_read_ids("SRR8073913.1 1", "SRR8073913.1 2"));
    assert!(compare_read_ids("SRR8073913.1/1", "SRR8073913.1/2"));

    // Case where IDs have common prefix but different SRA spot numbers (should NOT match)
    assert!(!compare_read_ids("SRR8073913.1.1", "SRR8073913.2.1"));
}

#[test]
fn test_complex_illumina_ids() {
    // NovaSeq / newer Illumina style
    assert!(compare_read_ids("A00123:123:ABC:1:1101:1234:5678 1:N:0:ATGC", "A00123:123:ABC:1:1101:1234:5678 2:N:0:ATGC"));
}
