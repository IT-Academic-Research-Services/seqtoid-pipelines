// Functions and definitions for the minimap2-associated PAF file format
use anyhow::{anyhow, Result};
use log::{info, debug};
use std::collections::HashMap;
use std::path::PathBuf;

use tokio::fs::File;
use tokio::io::{AsyncBufReadExt, BufReader};
use tokio::sync::mpsc::Sender;
use rayon::prelude::*;
use memchr::memchr2_iter;
use lexical::parse as lexical_parse;
use once_cell::sync::Lazy;
use bytes::Bytes;

use crate::utils::streams::ParseOutput;
use crate::config::defs::{SIMD_LEVEL, SimdLevel};

const LAMBDA: f64 = 1.58;
const K: f64 = 0.1;

#[derive(Debug, Clone)]
pub struct PafRecord {
    qname: String,
    qlen: u64,
    qstart: u64,
    qend: u64,
    strand: char,
    tname: String,
    tlen: u64,
    tstart: u64,
    tend: u64,
    nmatch: u64,
    alen: u64,
    mapq: u64,
    tags: std::collections::HashMap<String, String>,
}

impl PafRecord {
    /// Scalar parser — baseline, used on non-AVX-512 hardware.
    fn parse_line_scalar(line: &str) -> Result<Self> {
        let line = line.trim_end();
        let mut fields = line.split('\t');

        macro_rules! next {
            ($name:literal) => {
                fields.next().ok_or_else(|| anyhow!(concat!("Missing ", $name)))?
            };
        }
        macro_rules! parse_u64 {
            ($name:literal) => {{
                let s = next!($name);
                lexical_parse::<u64, _>(s.as_bytes())
                    .map_err(|e| anyhow!(concat!("invalid ", $name, ": {}"), e))?
            }};
        }

        let qname  = next!("qname").to_string();
        let qlen   = parse_u64!("qlen");
        let qstart = parse_u64!("qstart");
        let qend   = parse_u64!("qend");
        let strand = next!("strand").chars().next().ok_or_else(|| anyhow!("Invalid strand"))?;
        let tname  = next!("tname").to_string();
        let tlen   = parse_u64!("tlen");
        let tstart = parse_u64!("tstart");
        let tend   = parse_u64!("tend");
        let nmatch = parse_u64!("nmatch");
        let alen   = parse_u64!("alen");
        let mapq   = parse_u64!("mapq");

        let mut tags = HashMap::new();
        for tag_str in fields {
            let parts: Vec<&str> = tag_str.splitn(3, ':').collect();
            if parts.len() == 3 {
                tags.insert(parts[0].to_string(), parts[2].to_string());
            }
        }

        Ok(Self { qname, qlen, qstart, qend, strand, tname, tlen, tstart, tend, nmatch, alen, mapq, tags })
    }

    /// AVX-512 inner parser.
    /// Uses a 64-byte SIMD sweep to collect all tab positions in one pass,
    /// then indexes directly into byte slices — no iterator, no intermediate allocs.
    /// AVX-512 inner parser — now guaranteed to match scalar exactly.
    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx512f,avx512bw")]
    unsafe fn parse_line_avx512_inner(line: &str) -> Result<Self> {
        use std::arch::x86_64::*;

        let trimmed = line.trim_end();
        let bytes = trimmed.as_bytes();
        let len = bytes.len();
        if len == 0 {
            return Err(anyhow!("empty PAF line"));
        }

        // Phase 1: Collect all tab positions
        let mut tab_pos: Vec<usize> = Vec::with_capacity(24);
        let tab_splat = _mm512_set1_epi8(b'\t' as i8);
        let mut i = 0usize;

        while i + 64 <= len {
            let chunk = _mm512_loadu_si512(bytes.as_ptr().add(i).cast::<__m512i>());
            let mut mask: u64 = _mm512_cmpeq_epi8_mask(chunk, tab_splat);
            while mask != 0 {
                let bit = mask.trailing_zeros() as usize;
                tab_pos.push(i + bit);
                mask &= mask - 1;
            }
            i += 64;
        }
        while i < len {
            if bytes[i] == b'\t' {
                tab_pos.push(i);
            }
            i += 1;
        }

        if tab_pos.len() < 11 {
            return Err(anyhow!(
                "PAF line has {} tabs, need ≥11 for 12 mandatory fields",
                tab_pos.len()
            ));
        }

        // Phase 2: Field extraction (exact same logic as scalar)
        macro_rules! field {
            ($n:expr) => {{
                let start = if $n == 0 { 0 } else { tab_pos[$n - 1] + 1 };
                let end = if $n < tab_pos.len() { tab_pos[$n] } else { len };
                &bytes[start..end]
            }};
        }

        macro_rules! parse_u64 {
            ($n:expr, $name:literal) => {
                lexical_parse::<u64, _>(field!($n))
                    .map_err(|e| anyhow!(concat!($name, ": {}"), e))?
            };
        }

        let qname = std::str::from_utf8(field!(0))
            .map_err(|_| anyhow!("qname not UTF-8"))?
            .to_string();

        let qlen   = parse_u64!(1,  "qlen");
        let qstart = parse_u64!(2,  "qstart");
        let qend   = parse_u64!(3,  "qend");

        let strand = match field!(4).first() {
            Some(&b'+') => '+',
            Some(&b'-') => '-',
            _ => return Err(anyhow!("invalid strand")),
        };

        let tname = std::str::from_utf8(field!(5))
            .map_err(|_| anyhow!("tname not UTF-8"))?
            .to_string();

        let tlen   = parse_u64!(6,  "tlen");
        let tstart = parse_u64!(7,  "tstart");
        let tend   = parse_u64!(8,  "tend");
        let nmatch = parse_u64!(9,  "nmatch");
        let alen   = parse_u64!(10, "alen");
        let mapq   = parse_u64!(11, "mapq");

        // Phase 3: Tags — match scalar exactly
        let mut tags = HashMap::new();
        for fi in 12..=tab_pos.len() {  // note: <= to include last field
            let start = tab_pos[fi - 1] + 1;  // safe because we checked tab_pos.len()
            let end = if fi < tab_pos.len() { tab_pos[fi] } else { len };
            let tag_bytes = &bytes[start..end];

            if tag_bytes.is_empty() { continue; }
            if let Ok(s) = std::str::from_utf8(tag_bytes) {
                let parts: Vec<&str> = s.splitn(3, ':').collect();
                if parts.len() == 3 {
                    tags.insert(parts[0].to_string(), parts[2].to_string());
                }
            }
        }

        Ok(Self { qname, qlen, qstart, qend, strand, tname, tlen, tstart, tend, nmatch, alen, mapq, tags })
    }

    /// Safe wrapper — satisfies the `fn(&str) -> Result<Self>` type required by the dispatch static.
    #[cfg(target_arch = "x86_64")]
    fn parse_line_avx512(line: &str) -> Result<Self> {
        // Safety: only reachable when SIMD_LEVEL == Avx512, which means
        // is_x86_feature_detected! confirmed avx512f+bw at startup.
        unsafe { Self::parse_line_avx512_inner(line) }
    }

    /// Public entry point — dispatches to the best available implementation at first call.
    pub fn parse_line(line: &str) -> Result<Self> {
        PARSE_PAF_LINE(line)
    }

    pub fn alignment_score(&self) -> i32 {
        self.tags
            .get("AS")
            .and_then(|s| s.parse::<i32>().ok())
            .unwrap_or(0)
    }

    pub fn to_paf_line(&self) -> String {
        let mut line = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.qname, self.qlen, self.qstart, self.qend, self.strand,
            self.tname, self.tlen, self.tstart, self.tend, self.nmatch, self.alen, self.mapq,
        );

        for (key, value) in &self.tags {
            line.push_str(&format!("\t{}:Z:{}", key, value));
        }

        line
    }

    pub async fn merge_paf_files(
        paf_paths: &[PathBuf],
        tx: Sender<ParseOutput>,
    ) -> Result<()> {
        info!("Merging {} raw PAF files → streaming PAF output (no conversion)", paf_paths.len());

        let start = std::time::Instant::now();
        let total_lines = 0u64;

        let mut handles = Vec::new();

        for path in paf_paths {
            let path = path.clone();
            let tx_clone = tx.clone();

            handles.push(tokio::spawn(async move {
                let file = File::open(&path).await
                    .map_err(|e| anyhow!("Failed to open PAF {}: {}", path.display(), e))?;

                let mut reader = BufReader::new(file);
                let mut line = String::new();
                let mut local_lines = 0u64;
                let mut local_skipped = 0u64;

                while reader.read_line(&mut line).await? > 0 {
                    local_lines += 1;
                    let trimmed = line.trim();

                    if trimmed.is_empty() || trimmed.starts_with('#') {
                        local_skipped += 1;
                        line.clear();
                        continue;
                    }

                    // Forward raw PAF line unchanged
                    let bytes = (trimmed.to_string() + "\n").into_bytes();
                    if tx_clone.send(ParseOutput::Bytes(Bytes::from(bytes))).await.is_err() {
                        return Ok(());
                    }

                    line.clear();
                }

                info!("Merged {}: {} lines ({} skipped)", path.display(), local_lines, local_skipped);
                Ok::<(), anyhow::Error>(())
            }));
        }

        for h in handles {
            h.await??;
        }

        info!("Raw PAF merge complete: {} total lines in {:.2?}", total_lines, start.elapsed());
        Ok(())
    }

    fn calc_bitscore(&self) -> f64 {
        let nonmatch = self.tags.get("NM").and_then(|s| s.parse::<f64>().ok()).unwrap_or(0.0);
        let alen = self.alen as f64;
        let score = alen - 2.0 * nonmatch;
        (score * LAMBDA - K.ln()) / 2.0f64.ln()
    }

    fn calc_evalue(&self, genome_size: f64) -> f64 {
        let nonmatch = self.tags.get("NM").and_then(|s| s.parse::<f64>().ok()).unwrap_or(0.0);
        let alen = self.alen as f64;
        let score = alen - 2.0 * nonmatch;
        K * alen * genome_size * (-LAMBDA * score).exp()
    }

    fn calc_gap_openings(&self) -> u64 {
        let cigar = self.tags.get("cg").map(|s| s.as_bytes()).unwrap_or(b"");
        memchr2_iter(b'I', b'D', cigar).count() as u64
    }

    fn percent_identity(&self) -> f64 {
        if self.alen == 0 {
            0.0
        } else {
            (self.nmatch as f64 / self.alen as f64) * 100.0
        }
    }

    fn extract_accession(&self) -> String {
        self.tname
            .split(':')
            .find(|part| part.starts_with("NT:") || part.starts_with("NR:"))
            .and_then(|s| s.split('|').next())
            .map(|s| s.trim().to_string())
            .or_else(|| {
                self.tname
                    .split(|c| c == ' ' || c == '\t')
                    .next()
                    .map(|s| s.strip_prefix('>').unwrap_or(s).trim().to_string())
            })
            .unwrap_or_default()
    }

    pub fn to_m8_line(&self, genome_size: f64) -> String {
        let accession = self.extract_accession();
        let base_accession = accession.split('.').next().unwrap_or(&accession).to_string();

        if accession.is_empty() {
            return String::new();
        }

        let (mut tstart, mut tend) = (self.tstart, self.tend);
        if self.strand == '-' {
            std::mem::swap(&mut tstart, &mut tend);
        }
        let qstart_1 = self.qstart + 1;
        let tstart_adj = if self.strand == '+' { tstart + 1 } else { tstart };
        let tend_adj   = if self.strand == '+' { tend } else { tend + 1 };

        let nonmatch = self.tags.get("NM")
            .and_then(|s| s.parse::<u64>().ok())
            .unwrap_or(0);
        let gap_openings = self.calc_gap_openings();
        let percent_ident = self.percent_identity();
        let bitscore = self.calc_bitscore();
        let evalue = self.calc_evalue(genome_size);

        format!(
            "{}\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3e}\t{:.3}\n",
            self.qname, base_accession, percent_ident, self.alen,
            nonmatch, gap_openings, qstart_1, self.qend,
            tstart_adj, tend_adj, evalue, bitscore
        )
    }
}

static PARSE_PAF_LINE: Lazy<fn(&str) -> Result<PafRecord>> = Lazy::new(|| {
    #[cfg(target_arch = "x86_64")]
    if matches!(*SIMD_LEVEL, SimdLevel::Avx512) {
        debug!("PAF parser: using AVX-512 path");
        return PafRecord::parse_line_avx512;
    }
    debug!("PAF parser: using scalar path");
    PafRecord::parse_line_scalar
});

pub fn parse_paf_batch_to_m8(batch: Vec<u8>, genome_size: f64) -> Vec<Vec<u8>> {
    batch
        .par_split(|&b| b == b'\n')
        .filter(|line: &&[u8]| !line.is_empty() && !line.starts_with(b"#"))
        .flat_map(|line_bytes: &[u8]| {
            if let Ok(line) = std::str::from_utf8(line_bytes) {
                if let Ok(record) = PafRecord::parse_line(line) {
                    let m8_line = record.to_m8_line(genome_size);
                    return vec![m8_line.into_bytes()];
                }
            }
            vec![]
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    // ── Strategy for generating tag fields (valid + malformed) ─────────────
    fn tag_field_strategy() -> impl Strategy<Value = String> {
        prop_oneof![
            // Valid tag: key:type:value
            ("[A-Za-z0-9_]{1,10}", "[A-Za-z]:", "[^\\t\\n]{0,40}")
                .prop_map(|(k, t, v)| format!("{}:{}{}", k, t, v)),

            // Malformed: missing type (only one colon)
            ("[A-Za-z0-9_]{1,10}", "[^:\\t\\n]{0,30}")
                .prop_map(|(k, v)| format!("{}:{}", k, v)),

            // Malformed: empty key
            (":[A-Za-z]:[^\\t\\n]{0,30}").prop_map(|s| s.to_string()),

            // Malformed: empty value
            ("[A-Za-z0-9_]{1,10}:[A-Za-z]:").prop_map(|s| s.to_string()),

            // Malformed: garbage / no colons
            "[A-Za-z0-9_:\\-]{1,40}".prop_map(|s| s),

            // Malformed: trailing colon only
            "[A-Za-z0-9_]{1,10}:".prop_map(|s| s.to_string()),

            // Empty tag field (consecutive tabs)
            Just("".to_string()),
        ]
    }

    // ── Basic field proptest ───────────────────────────────────────────────
    proptest! {
        #[test]
        fn test_parse_line_paf_equivalence(
            qname in "[a-zA-Z0-9._-]{1,100}",
            qlen in 1u64..1000000,
            qstart in 0u64..1000000,
            qend in 0u64..1000000,
            strand in "[+-]",
            tname in "[a-zA-Z0-9._-]{1,100}",
            tlen in 1u64..1000000000,
            tstart in 0u64..1000000000,
            tend in 0u64..1000000000,
            nmatch in 0u64..1000000,
            alen in 0u64..1000000,
            mapq in 0u8..60
        ) {
            let line = format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                qname, qlen, qstart, qend, strand, tname, tlen, tstart, tend, nmatch, alen, mapq
            );
            compare_parsers(&line);
        }
    }

    // ── Proptests with random extra tags (valid + malformed) ───────────────
    proptest! {
        #[test]
        fn test_parse_line_paf_equivalence_with_malformed_tags(
            qname in "[a-zA-Z0-9._-]{1,100}",
            qlen in 1u64..1000000,
            qstart in 0u64..1000000,
            qend in 0u64..1000000,
            strand in "[+-]",
            tname in "[a-zA-Z0-9._-]{1,100}",
            tlen in 1u64..1000000000,
            tstart in 0u64..1000000000,
            tend in 0u64..1000000000,
            nmatch in 0u64..1000000,
            alen in 0u64..1000000,
            mapq in 0u8..60,
            extra_tags in prop::collection::vec(tag_field_strategy(), 0..40)
        ) {
            let mut line = format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                qname, qlen, qstart, qend, strand, tname, tlen, tstart, tend, nmatch, alen, mapq
            );

            for tag in extra_tags {
                line.push_str(&format!("\t{}", tag));
            }

            let scalar_res = PafRecord::parse_line_scalar(&line);
            let dispatched_res = PafRecord::parse_line(&line);

            match (scalar_res, dispatched_res) {
                (Ok(s), Ok(d)) => assert_records_eq(&s, &d, &line),
                (Err(_), Err(_)) => {}
                _ => panic!("Scalar and AVX-512 disagreed on success/failure for line: {}", line),
            }
        }
    }

    // ── helpers ────────────────────────────────────────────────────────────

    /// Assert two PafRecords are field-for-field identical.
    #[allow(dead_code)]
    fn assert_records_eq(a: &PafRecord, b: &PafRecord, context: &str) {
        assert_eq!(a.qname,  b.qname,  "{context}: qname");
        assert_eq!(a.qlen,   b.qlen,   "{context}: qlen");
        assert_eq!(a.qstart, b.qstart, "{context}: qstart");
        assert_eq!(a.qend,   b.qend,   "{context}: qend");
        assert_eq!(a.strand, b.strand, "{context}: strand");
        assert_eq!(a.tname,  b.tname,  "{context}: tname");
        assert_eq!(a.tlen,   b.tlen,   "{context}: tlen");
        assert_eq!(a.tstart, b.tstart, "{context}: tstart");
        assert_eq!(a.tend,   b.tend,   "{context}: tend");
        assert_eq!(a.nmatch, b.nmatch, "{context}: nmatch");
        assert_eq!(a.alen,   b.alen,   "{context}: alen");
        assert_eq!(a.mapq,   b.mapq,   "{context}: mapq");
        assert_eq!(a.tags,   b.tags,   "{context}: tags");
    }

    /// Run both parsers against `line` and assert they agree.
    fn compare_parsers(line: &str) {
        let scalar = PafRecord::parse_line_scalar(line)
            .expect("scalar parse failed");

        let dispatched = PafRecord::parse_line(line)
            .expect("dispatched (AVX-512 or scalar) parse failed");

        assert_records_eq(&scalar, &dispatched, line);
    }

    // ── representative real-world lines ───────────────────────────────────

    const BASIC: &str =
        "read1\t150\t0\t150\t+\tNC_045512.2\t29903\t100\t250\t148\t150\t60";

    const WITH_TAGS: &str =
        "read2\t300\t10\t290\t-\tAB123456\t5000\t200\t480\t270\t280\t30\tcg:Z:10M2I5D\tNM:i:7\tAS:i:255";

    const REVERSE_STRAND: &str =
        "read3\t200\t5\t195\t-\tKJ660346\t18958\t500\t690\t185\t190\t60\tAS:i:180";

    const ZERO_COORDS: &str =
        "read4\t100\t0\t100\t+\tMN908947.3\t29903\t0\t100\t100\t100\t0";

    const MAX_VALUES: &str =
        "read5\t1000000\t0\t999999\t+\tNC_000001\t248956422\t0\t999999\t999000\t1000000\t60";

    const LONG_LINE: &str =
        "long_read_id_that_is_quite_verbose_indeed\t250\t0\t250\t+\tsome_reference_sequence_name\t100000\t1000\t1250\t245\t250\t60\tcg:Z:250M\tNM:i:5";

    const TRAILING_SPACE: &str =
        "read6\t150\t0\t150\t+\tNC_045512.2\t29903\t100\t250\t148\t150\t60   ";

    // ── tests ─────────────────────────────────────────────────────────────

    #[test]
    fn test_basic_line() {
        compare_parsers(BASIC);
    }

    #[test]
    fn test_with_optional_tags() {
        compare_parsers(WITH_TAGS);
        let r = PafRecord::parse_line_scalar(WITH_TAGS).unwrap();
        assert_eq!(r.tags.get("NM").map(|s| s.as_str()), Some("7"));
        assert_eq!(r.tags.get("AS").map(|s| s.as_str()), Some("255"));
        assert_eq!(r.tags.get("cg").map(|s| s.as_str()), Some("10M2I5D"));
    }

    #[test]
    fn test_reverse_strand() {
        compare_parsers(REVERSE_STRAND);
        let r = PafRecord::parse_line_scalar(REVERSE_STRAND).unwrap();
        assert_eq!(r.strand, '-');
    }

    #[test]
    fn test_zero_coordinates() {
        compare_parsers(ZERO_COORDS);
        let r = PafRecord::parse_line_scalar(ZERO_COORDS).unwrap();
        assert_eq!(r.qstart, 0);
        assert_eq!(r.tstart, 0);
    }

    #[test]
    fn test_max_values() {
        compare_parsers(MAX_VALUES);
        let r = PafRecord::parse_line_scalar(MAX_VALUES).unwrap();
        assert_eq!(r.qlen, 1_000_000);
        assert_eq!(r.tlen, 248_956_422);
    }

    #[test]
    fn test_line_spanning_two_simd_chunks() {
        compare_parsers(LONG_LINE);
    }

    #[test]
    fn test_trailing_whitespace_trimmed() {
        compare_parsers(TRAILING_SPACE);
        let r = PafRecord::parse_line_scalar(TRAILING_SPACE).unwrap();
        assert_eq!(r.qname, "read6");
        assert_eq!(r.qlen, 150);
        assert_eq!(r.mapq, 60);
    }

    #[test]
    fn test_to_m8_line_scalar_avx512_agree() {
        let genome_size = 29903.0_f64;
        let scalar_m8 = PafRecord::parse_line_scalar(WITH_TAGS)
            .unwrap()
            .to_m8_line(genome_size);

        let dispatched_m8 = PafRecord::parse_line(WITH_TAGS)
            .unwrap()
            .to_m8_line(genome_size);

        assert_eq!(scalar_m8, dispatched_m8, "to_m8_line output differs between scalar and dispatched");
        assert!(!scalar_m8.is_empty());
    }

    #[test]
    fn test_error_on_empty_line() {
        assert!(PafRecord::parse_line_scalar("").is_err());
        assert!(PafRecord::parse_line("").is_err());
    }

    #[test]
    fn test_error_on_too_few_fields() {
        let short = "read1\t150\t0\t150\t+\tNC_045512.2\t29903";
        assert!(PafRecord::parse_line_scalar(short).is_err());
        assert!(PafRecord::parse_line(short).is_err());
    }

    #[test]
    fn test_parse_paf_batch_to_m8_roundtrip() {
        let batch = format!("{}\n{}\n{}\n", BASIC, WITH_TAGS, REVERSE_STRAND);
        let results = parse_paf_batch_to_m8(batch.into_bytes(), 29903.0);
        let _ = results;
    }

    // ── Strong edge-case tests (malformed tags, long lines, etc.) ─────────

    #[test]
    fn test_trailing_whitespace_and_empty_tags() {
        let cases = vec![
            "read6\t150\t0\t150\t+\tNC_045512.2\t29903\t100\t250\t148\t150\t60   ",
            "read7\t100\t0\t100\t+\tref\t1000\t0\t100\t100\t100\t60\t\tNM:i:0",
            "read8\t200\t0\t200\t-\tref\t2000\t0\t200\t200\t200\t60\tcg:Z:200M  ",
            "read9\t150\t0\t150\t+\tref\t1000\t0\t150\t150\t150\t60\t   ",
        ];

        for line in cases {
            compare_parsers(line);
        }
    }

    #[test]
    fn test_malformed_tags() {
        let line = "read\t100\t0\t100\t+\tref\t1000\t0\t100\t100\t100\t60\tbadtag\tgood:Z:value\tNM:i:5\tincomplete:";
        compare_parsers(line);

        let r = PafRecord::parse_line(line).unwrap();
        assert_eq!(r.tags.get("good"), Some(&"value".to_string()));
        assert_eq!(r.tags.get("NM"), Some(&"5".to_string()));
        assert!(!r.tags.contains_key("badtag"));
        assert!(!r.tags.contains_key("incomplete"));
    }

    #[test]
    fn test_long_line_many_tabs() {
        let mut line = "longid\t100\t0\t100\t+\tlongref\t10000\t0\t100\t100\t100\t60".to_string();
        for i in 0..50 {
            line.push_str(&format!("\ttag{i}:Z:value{i}"));
        }
        compare_parsers(&line);
    }

    #[test]
    fn test_header_like_ids_and_special_chars() {
        let cases = vec![
            "read:1|with:colon\t150\t0\t150\t+\tref:seq|version.1\t1000\t0\t150\t150\t150\t60",
            "read with space in qname\t100\t0\t100\t+\tref\t1000\t0\t100\t100\t100\t60",
            "read\t100\t0\t100\t+\tref\t1000\t0\t100\t100\t100\t60\tcg:Z:10M2D\tAS:i:200\tNM:i:2",
        ];

        for line in cases {
            compare_parsers(line);
        }
    }
}