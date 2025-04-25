use rand::rngs::ThreadRng;
use rand::seq::IndexedRandom;
use rand::rng;

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
        *DNA::all()
            .choose(rng)
            .expect("DNA::all is never empty")
    }
}

fn phred33(score: u8) -> u8 {
    score + 33
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_random_nucleotide() {
        let dna = DNA::random();
        assert!(matches!(
            dna,
            DNA::A | DNA::C | DNA::G | DNA::T
        ));
    }

    #[test]
    fn test_random_sequence() {
        let seq = DNA::random_sequence(10);
        assert_eq!(seq.len(), 10);
        assert!(seq.chars().all(|c| "ACGT".contains(c)));
    }
}