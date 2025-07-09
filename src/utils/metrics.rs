/// Calculation of metrics from result files.


// MUMmer
// For the show-coords executable from MUMmer suite
#[derive(Debug)]
pub struct AlignmentMetrics {
    pub total_aligned_length: u64,
    pub largest_alignment: u64,
    pub genome_fraction: f64,
    pub na50: u64,
    pub nga50: u64,
    pub na75: u64,
    pub nga75: u64,
    pub la50: u32,
    pub lga50: u32,
    pub la75: u32,
    pub lga75: u32,
    pub mismatches_per_100kbp: f64,
    pub indels_per_100kbp: f64,
    pub misassemblies: u32,
    pub unaligned_length: u64,
    pub duplication_ratio: f64,
}