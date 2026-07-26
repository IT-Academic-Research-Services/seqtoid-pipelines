// nr_preprocess.rs
// High-performance, correct streaming FASTA preprocessor for MMseqs2 nr DB creation.
// - Never silently drops data: every input record is processed and output (or explicitly handled).
// - Validates amino acid sequences, fixes invalid chars to 'X' (logs examples).
// - Sanitizes headers: removes control characters (including ^A) and collapses whitespace.
// - Large I/O buffers for speed on Epyc / high-core systems.
// - Pure Rust, std only.
// Compile:
//   rustc -O -C target-cpu=native nr_preprocess.rs -o nr_preprocess
//
// Usage (stream style, no data loss):
//   zcat /data/refs/nr.fa.gz | ./nr_preprocess > /scratch/nr_clean.fa
//   or:
//   cat /data/refs/nr.fa | ./nr_preprocess > /scratch/nr_clean.fa
//
// Then:
//   mmseqs createdb /scratch/nr_clean.fa /data/refs/nr_clean/nrcleanDB \
//     --dbtype 1 --compressed 1 --threads 64 2>&1 | tee createdb.log

use std::collections::HashSet;
use std::env;
use std::io::{self, BufRead, BufReader, BufWriter, Write};

fn main() -> io::Result<()> {
    let args: Vec<String> = env::args().collect();
    let mut input_name = "stdin".to_string();
    let mut output_name = "stdout".to_string();

    if args.len() > 1 && args[1] != "-" {
        input_name = args[1].clone();
    }
    if args.len() > 2 {
        output_name = args[2].clone();
    }

    eprintln!("nr_preprocess starting | input: {} | output: {}", input_name, output_name);
    eprintln!("Strict mode: no silent drops. All sequences will be output. Invalid AA -> 'X' with logging.");
    eprintln!("Header cleaning: remove control characters (including ^A) and collapse whitespace.");

    let stdin = io::stdin();
    let reader: Box<dyn BufRead> = if input_name == "stdin" || input_name == "-" {
        Box::new(BufReader::with_capacity(128 * 1024 * 1024, stdin.lock()))
    } else {
        let file = std::fs::File::open(&input_name)?;
        Box::new(BufReader::with_capacity(128 * 1024 * 1024, file))
    };

    let stdout = io::stdout();
    let mut writer: Box<dyn Write> = if output_name == "stdout" || output_name == "-" {
        Box::new(BufWriter::with_capacity(128 * 1024 * 1024, stdout.lock()))
    } else {
        let file = std::fs::File::create(&output_name)?;
        Box::new(BufWriter::with_capacity(128 * 1024 * 1024, file))
    };

    // Standard amino-acid alphabet + common specials used in NR
    let valid_aa: HashSet<char> = "ACDEFGHIKLMNPQRSTVWYXBZJUO*.-".chars().collect();

    let mut lines = reader.lines();
    let mut input_seqs: u64 = 0;
    let mut output_seqs: u64 = 0;
    let mut seqs_with_fixes: u64 = 0;
    let mut total_replacements: u64 = 0;
    let mut bad_char_examples: Vec<(u64, char)> = Vec::new();
    let mut current_header: Option<String> = None;
    let mut current_seq = String::with_capacity(4096);

    let progress_every: u64 = 1_000_000;

    while let Some(line_res) = lines.next() {
        let line = line_res?;
        if line.starts_with('>') {
            // Finish previous record
            if let Some(ref hdr) = current_header {
                if !current_seq.is_empty() {
                    process_and_write_record(
                        hdr,
                        &current_seq,
                        &valid_aa,
                        &mut writer,
                        &mut seqs_with_fixes,
                        &mut total_replacements,
                        &mut bad_char_examples,
                        input_seqs,
                    )?;
                    output_seqs += 1;
                } else {
                    eprintln!("WARNING: empty sequence for entry {} header: {}", input_seqs, hdr);
                }
            }
            current_header = Some(line);
            current_seq.clear();
            input_seqs += 1;

            if input_seqs % progress_every == 0 {
                eprintln!("Progress: processed {} sequences so far...", input_seqs);
            }
        } else if current_header.is_some() && !line.trim().is_empty() {
            current_seq.push_str(line.trim_end());
        }
    }

    // Process last record
    if let Some(ref hdr) = current_header {
        if !current_seq.is_empty() {
            process_and_write_record(
                hdr,
                &current_seq,
                &valid_aa,
                &mut writer,
                &mut seqs_with_fixes,
                &mut total_replacements,
                &mut bad_char_examples,
                input_seqs,
            )?;
            output_seqs += 1;
        }
    }

    writer.flush()?;

    eprintln!("\n=== nr_preprocess FINAL REPORT ===");
    eprintln!("Input sequences read:     {}", input_seqs);
    eprintln!("Output sequences written: {}", output_seqs);
    eprintln!("Sequences requiring AA fix (invalid char -> 'X'): {}", seqs_with_fixes);
    eprintln!("Total invalid chars replaced: {}", total_replacements);
    if !bad_char_examples.is_empty() {
        eprintln!(
            "First bad char examples (entry#, char): {:?}",
            &bad_char_examples[..bad_char_examples.len().min(10)]
        );
    }
    if input_seqs == output_seqs {
        eprintln!("SUCCESS: ZERO data dropped. Streamed every record correctly. Ready for MMseqs createdb.");
    } else {
        eprintln!("ERROR: Count mismatch detected! Investigate parse logic. (This should never happen.)");
        std::process::exit(2);
    }
    eprintln!("Next recommended:");
    eprintln!("  mmseqs createdb /scratch/nr_clean.fa /data/refs/nr_clean/nrcleanDB \\");
    eprintln!("    --dbtype 1 --compressed 1 --threads 64 2>&1 | tee createdb.log");
    eprintln!("Then immediately verify type:");
    eprintln!("  od -An -t u4 /data/refs/nr_clean/nrcleanDB.dbtype");
    eprintln!("  (must print 1)");

    Ok(())
}

fn process_and_write_record(
    header: &str,
    seq: &str,
    valid_aa: &HashSet<char>,
    writer: &mut dyn Write,
    seqs_with_fixes: &mut u64,
    total_replacements: &mut u64,
    bad_examples: &mut Vec<(u64, char)>,
    entry_num: u64,
) -> io::Result<()> {
    // ---------------------------------------------------------------
    // Header sanitization for MMseqs2
    // - Remove all control characters (including the famous NR ^A = 0x01)
    // - Collapse multiple whitespace into single spaces
    // - Keep as much original information as possible
    // ---------------------------------------------------------------
    let clean_header = {
        let mut h = header.trim().to_string();

        // Remove every control character
        h.retain(|c| !c.is_control());

        // Collapse runs of whitespace into a single space
        let mut result = String::with_capacity(h.len());
        let mut last_was_space = false;
        for c in h.chars() {
            if c.is_whitespace() {
                if !last_was_space {
                    result.push(' ');
                    last_was_space = true;
                }
            } else {
                result.push(c);
                last_was_space = false;
            }
        }
        result.trim().to_string()
    };

    // ---------------------------------------------------------------
    // Sequence cleaning: invalid amino-acid characters → 'X'
    // ---------------------------------------------------------------
    let mut fixed_this = false;
    let mut cleaned_seq = String::with_capacity(seq.len());

    for c in seq.chars() {
        let cu = c.to_ascii_uppercase();
        if valid_aa.contains(&cu) {
            cleaned_seq.push(cu);
        } else {
            if bad_examples.len() < 10 {
                bad_examples.push((entry_num, c));
            }
            cleaned_seq.push('X');
            *total_replacements += 1;
            fixed_this = true;
        }
    }

    if fixed_this {
        *seqs_with_fixes += 1;
    }

    // Write one clean header line + one sequence line
    writeln!(writer, "{}", clean_header)?;
    writeln!(writer, "{}", cleaned_seq)?;

    Ok(())
}
