//! Sequence utilities.

use rand::rng;
use rand::rngs::ThreadRng;
use rand::seq::IndexedRandom;
use rand_distr::{Distribution, Normal};

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DNA {
    A,
    C,
    G,
    T,
}

impl DNA {
    /// Convert nucleotide to its character representation.
    pub fn to_char(&self) -> char {
        match self {
            DNA::A => 'A',
            DNA::C => 'C',
            DNA::G => 'G',
            DNA::T => 'T',
        }
    }

    /// Get all possible nucleotides as a static slice.
    pub fn all() -> &'static [DNA] {
        &[DNA::A, DNA::C, DNA::G, DNA::T]
    }

    /// Generate a random nucleotide using the thread-local RNG. words.choose(&mut rand::rng()).unwrap()
    #[allow(dead_code)]
    pub fn random() -> DNA {
        let mut rng = rng();

        *DNA::all()
            .choose(&mut rng)
            .expect("Nucleotide::all is never empty")
    }

    /// Generate a random sequence of nucleotides of the given length.
    pub fn random_sequence(length: usize) -> String {
        let mut rng = rng();
        (0..length)
            .map(|_| DNA::random_with_rng(&mut rng).to_char())
            .collect()
    }

    /// Helper method to generate a random nucleotide with a provided RNG.
    fn random_with_rng(rng: &mut ThreadRng) -> DNA {
        *DNA::all().choose(rng).expect("DNA::all is never empty")
    }
}

pub mod valid_bases {
    /// Valid DNA nucleotide bases in ASCII: A, C, G, T, a, c, g, t
    pub const DNA: &[u8] = &[65, 67, 71, 84, 97, 99, 103, 116];

    /// Valid DNA nucleotide bases including 'N' in ASCII: A, C, G, T, N, a, c, g, t, n
    pub const DNA_WITH_N: &[u8] = &[65, 67, 71, 84, 78, 97, 99, 103, 116, 110];
}


/// Error returned when a DNA sequence cannot be translated with the
/// RGI/Biopython-compatible NCBI translation-table-11 rules.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TranslationError(pub String);

impl std::fmt::Display for TranslationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for TranslationError {}

const TABLE11_START_CODONS: [&[u8; 3]; 7] = [
    b"TTG", b"CTG", b"ATT", b"ATC", b"ATA", b"ATG", b"GTG",
];
const TABLE11_STOP_CODONS: [&[u8; 3]; 3] = [b"TAA", b"TAG", b"TGA"];

fn table11_unambiguous_codon(codon: &[u8; 3]) -> Option<char> {
    Some(match codon {
        b"TTT" | b"TTC" => 'F',
        b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG" => 'L',
        b"ATT" | b"ATC" | b"ATA" => 'I',
        b"ATG" => 'M',
        b"GTT" | b"GTC" | b"GTA" | b"GTG" => 'V',
        b"TCT" | b"TCC" | b"TCA" | b"TCG" | b"AGT" | b"AGC" => 'S',
        b"CCT" | b"CCC" | b"CCA" | b"CCG" => 'P',
        b"ACT" | b"ACC" | b"ACA" | b"ACG" => 'T',
        b"GCT" | b"GCC" | b"GCA" | b"GCG" => 'A',
        b"TAT" | b"TAC" => 'Y',
        b"TAA" | b"TAG" | b"TGA" => '*',
        b"CAT" | b"CAC" => 'H',
        b"CAA" | b"CAG" => 'Q',
        b"AAT" | b"AAC" => 'N',
        b"AAA" | b"AAG" => 'K',
        b"GAT" | b"GAC" => 'D',
        b"GAA" | b"GAG" => 'E',
        b"TGT" | b"TGC" => 'C',
        b"TGG" => 'W',
        b"CGT" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG" => 'R',
        b"GGT" | b"GGC" | b"GGA" | b"GGG" => 'G',
        _ => return None,
    })
}

fn iupac_bases(base: u8) -> Option<&'static [u8]> {
    Some(match base.to_ascii_uppercase() {
        b'A' => b"A",
        b'C' => b"C",
        b'G' => b"G",
        b'T' | b'U' => b"T",
        b'R' => b"AG",
        b'Y' => b"CT",
        b'S' => b"GC",
        b'W' => b"AT",
        b'K' => b"GT",
        b'M' => b"AC",
        b'B' => b"CGT",
        b'D' => b"AGT",
        b'H' => b"ACT",
        b'V' => b"ACG",
        b'N' => b"ACGT",
        _ => return None,
    })
}

/// Translate one possibly-ambiguous codon using NCBI table 11.
///
/// Match Biopython's ambiguous-codon translation: return the unique amino
/// acid when possible, B for D/N, Z for E/Q, J for I/L, and X when no
/// narrower extended-protein symbol represents all possible translations.
fn table11_codon(codon: &[u8; 3]) -> Result<char, TranslationError> {
    let b0 = iupac_bases(codon[0]).ok_or_else(|| {
        TranslationError(format!("invalid nucleotide '{}'", codon[0] as char))
    })?;
    let b1 = iupac_bases(codon[1]).ok_or_else(|| {
        TranslationError(format!("invalid nucleotide '{}'", codon[1] as char))
    })?;
    let b2 = iupac_bases(codon[2]).ok_or_else(|| {
        TranslationError(format!("invalid nucleotide '{}'", codon[2] as char))
    })?;

    let mut translated: Option<char> = None;
    for &a in b0 {
        for &b in b1 {
            for &c in b2 {
                let concrete = [a, b, c];
                let aa = table11_unambiguous_codon(&concrete)
                    .expect("expanded IUPAC codon must be A/C/G/T");
                match translated {
                    None => translated = Some(aa),
                    Some(previous) if previous == aa => {}
                    Some(previous) => {
                        // Match Biopython's extended ambiguous-protein alphabet.
                        // B = D/N, Z = E/Q, J = I/L. Any other mixture,
                        // including STOP mixed with an amino acid, is X.
                        translated = Some(match (previous, aa) {
                            ('B', 'D' | 'N') | ('D' | 'N', 'B') => 'B',
                            ('D', 'N') | ('N', 'D') => 'B',
                            ('Z', 'E' | 'Q') | ('E' | 'Q', 'Z') => 'Z',
                            ('E', 'Q') | ('Q', 'E') => 'Z',
                            ('J', 'I' | 'L') | ('I' | 'L', 'J') => 'J',
                            ('I', 'L') | ('L', 'I') => 'J',
                            _ => 'X',
                        });

                        if translated == Some('X') {
                            return Ok('X');
                        }
                    }
                }
            }
        }
    }
    Ok(translated.expect("IUPAC expansion is never empty"))
}

fn translate_table11(seq: &[u8], cds: bool) -> Result<String, TranslationError> {
    if seq.len() % 3 != 0 {
        return Err(TranslationError(format!(
            "sequence length {} is not a multiple of three",
            seq.len()
        )));
    }
    if seq.is_empty() {
        return Ok(String::new());
    }

    let upper: Vec<u8> = seq.iter().map(|b| b.to_ascii_uppercase()).collect();

    if cds {
        let first: [u8; 3] = upper[0..3].try_into().unwrap();
        if !TABLE11_START_CODONS.iter().any(|start| **start == first) {
            return Err(TranslationError(format!(
                "first codon '{}' is not a table-11 start codon",
                String::from_utf8_lossy(&first)
            )));
        }

        let last: [u8; 3] = upper[upper.len() - 3..].try_into().unwrap();
        if !TABLE11_STOP_CODONS.iter().any(|stop| **stop == last) {
            return Err(TranslationError(format!(
                "final codon '{}' is not a table-11 stop codon",
                String::from_utf8_lossy(&last)
            )));
        }
    }

    let codon_count = upper.len() / 3;
    let mut protein = String::with_capacity(codon_count);
    for (index, chunk) in upper.chunks_exact(3).enumerate() {
        let codon: [u8; 3] = chunk.try_into().unwrap();

        // Biopython cds=True translates any valid table-11 initiation codon
        // as methionine and excludes the required terminal stop from output.
        if cds && index == 0 {
            protein.push('M');
            continue;
        }
        if cds && index + 1 == codon_count {
            break;
        }

        let aa = table11_codon(&codon)?;
        if cds && aa == '*' {
            return Err(TranslationError(format!(
                "internal stop codon at codon {}",
                index + 1
            )));
        }
        protein.push(aa);
    }

    Ok(protein)
}

/// Translate a KMA consensus exactly as RGI's BWT.py does:
/// pad an incomplete terminal codon with N, try NCBI table 11 with CDS
/// validation, fall back to ordinary table-11 translation on CDS failure,
/// then remove the amino acid produced by the artificial padded codon.
pub fn translate_rgi_kma_consensus(dna: &str) -> Result<String, TranslationError> {
    let mut input = dna.as_bytes().to_vec();
    let trailing_bases = input.len() % 3;

    if trailing_bases == 1 {
        input.extend_from_slice(b"NN");
    } else if trailing_bases == 2 {
        input.push(b'N');
    }

    let mut protein = match translate_table11(&input, true) {
        Ok(protein) => protein,
        Err(_) => translate_table11(&input, false)?,
    };

    if trailing_bases != 0 {
        protein.pop();
    }

    Ok(protein)
}

fn phred33(score: u8) -> u8 {
    score + 33
}

fn normal_phred_qual(mean: f32, stdev: f32) -> u8 {
    let mut raw_phred = -1.0;

    let normal = Normal::new(mean, stdev).unwrap();

    while raw_phred < 0.0 || raw_phred > 40.0 {
        raw_phred = normal.sample(&mut rand::rng());
    }

    phred33(raw_phred as u8)
}

pub fn normal_phred_qual_string(length: usize, mean: f32, stdev: f32) -> String {
    let mut quals = String::new();

    for _i in 0..length {
        quals.push(normal_phred_qual(mean, stdev) as char);
    }

    quals
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_table11_cds_alternative_start_is_methionine() {
        assert_eq!(translate_table11(b"GTGAAATAA", true).unwrap(), "MK");
    }

    #[test]
    fn test_table11_cds_failure_falls_back_for_rgi() {
        // AAA is not a table-11 initiation codon. RGI catches the CDS
        // translation failure and retranslates codon-by-codon.
        assert_eq!(translate_rgi_kma_consensus("AAAGGG").unwrap(), "KG");
    }

    #[test]
    fn test_rgi_translation_removes_padded_terminal_codon() {
        assert_eq!(translate_rgi_kma_consensus("ATGAA").unwrap(), "M");
    }

    #[test]
    fn test_ambiguous_codon_translation() {
        assert_eq!(table11_codon(b"GCN").unwrap(), 'A');
        assert_eq!(table11_codon(b"RAY").unwrap(), 'B'); // AAY/AAC? RAY -> D/N
        assert_eq!(table11_codon(b"SAR").unwrap(), 'Z'); // CAG/GAG -> Q/E
        assert_eq!(table11_codon(b"WTA").unwrap(), 'J'); // ATA/TTA -> I/L
        assert_eq!(table11_codon(b"TAN").unwrap(), 'X'); // Y or STOP
        assert_eq!(table11_codon(b"NNN").unwrap(), 'X');
    }

    #[test]
    fn test_random_nucleotide() {
        let dna = DNA::random();
        assert!(matches!(dna, DNA::A | DNA::C | DNA::G | DNA::T));
    }

    #[test]
    fn test_random_sequence() {
        let seq = DNA::random_sequence(10);
        assert_eq!(seq.len(), 10);
        assert!(seq.chars().all(|c| "ACGT".contains(c)));
    }
}
