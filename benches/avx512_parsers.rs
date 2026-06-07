use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use seqtoid_pipelines::utils::fastx::{fastx_generator, SequenceRecord, parse_header};
use std::io::Cursor;
use needletail::parser::FastqReader;
use needletail::FastxReader;

use seqtoid_pipelines::utils::blast::M8Record;
use seqtoid_pipelines::utils::paf::PafRecord;


// ── Existing micro benchmarks (kept for reference) ───────────────────────

fn bench_parse_header(c: &mut Criterion) {
    let long_header: &[u8] = b"very_long_read_id_that_spans_multiple_64_byte_chunks_xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx some description here";

    let mut group = c.benchmark_group("parse_header");
    group.throughput(Throughput::Bytes(long_header.len() as u64));

    group.bench_function("scalar", |b| {
        b.iter(|| parse_header(black_box(long_header), '>'))
    });

    group.bench_function("dispatched", |b| {
        b.iter(|| seqtoid_pipelines::utils::fastx::parse_header(black_box(long_header), '>'))
    });
    group.finish();
}

fn bench_paf_line(c: &mut Criterion) {
    let line = "read1\t150\t0\t150\t+\tNC_045512.2\t29903\t100\t250\t148\t150\t60\tcg:Z:150M\tNM:i:0\tAS:i:300";

    let mut group = c.benchmark_group("paf_single_line");
    group.throughput(Throughput::Bytes(line.len() as u64));

    group.bench_function("scalar", |b| {
        b.iter(|| PafRecord::parse_line_scalar(black_box(line)))
    });

    group.bench_function("dispatched", |b| {
        b.iter(|| PafRecord::parse_line(black_box(line)))
    });
    group.finish();
}

fn bench_m8_nt_line(c: &mut Criterion) {
    let line = "read1\tNC_045512.2\t99.333\t150\t1\t0\t1\t150\t100\t249\t1.23e-75\t285.000\t150\t29903\tAS:i:300\tNM:i:1";

    let mut group = c.benchmark_group("m8_nt_single_line");
    group.throughput(Throughput::Bytes(line.len() as u64));

    group.bench_function("scalar", |b| {
        b.iter(|| M8Record::parse_line_nt_scalar(black_box(line)))
    });

    group.bench_function("dispatched", |b| {
        b.iter(|| M8Record::parse_line_nt(black_box(line)))
    });
    group.finish();
}

fn bench_m8_nr_line(c: &mut Criterion) {
    let line = "read1\tQIK02963.1\t87.500\t96\t12\t0\t1\t96\t1\t96\t1.45e-38\t152.000\tAS:i:180\tNM:i:3";

    let mut group = c.benchmark_group("m8_nr_single_line");
    group.throughput(Throughput::Bytes(line.len() as u64));

    group.bench_function("scalar", |b| {
        b.iter(|| M8Record::parse_line_nr_scalar(black_box(line)))
    });

    group.bench_function("dispatched", |b| {
        b.iter(|| M8Record::parse_line_nr(black_box(line)))
    });
    group.finish();
}

// ── NEW: Large-Batch Streaming Benchmarks ────────────────────────────────

/// Generate a realistic large batch of NT M8 lines in memory (no disk I/O)
fn generate_m8_nt_batch(count: usize) -> String {
    let mut batch = String::with_capacity(count * 120);
    for i in 0..count {
        batch.push_str(&format!(
            "read{}\tNC_045512.2\t99.333\t150\t1\t0\t1\t150\t100\t249\t1.23e-75\t285.000\t150\t29903\n",
            i
        ));
    }
    batch
}

/// Generate a realistic large batch of PAF lines in memory
fn generate_paf_batch(count: usize) -> String {
    let mut batch = String::with_capacity(count * 110);
    for i in 0..count {
        batch.push_str(&format!(
            "read{}\t150\t0\t150\t+\tNC_045512.2\t29903\t100\t250\t148\t150\t60\tcg:Z:150M\tNM:i:0\n",
            i
        ));
    }
    batch
}

fn bench_m8_nt_large_batch_streaming(c: &mut Criterion) {
    // 500k lines ≈ realistic mNGS BLAST output size
    let batch = generate_m8_nt_batch(500_000);
    let bytes = batch.as_bytes();

    let mut group = c.benchmark_group("m8_nt_large_batch_streaming");
    group.throughput(Throughput::Bytes(bytes.len() as u64));
    group.sample_size(10); // fewer samples because this is heavy

    group.bench_function("scalar", |b| {
        b.iter(|| {
            for line in bytes.split(|&c| c == b'\n') {
                if line.is_empty() { continue; }
                let _ = black_box(M8Record::parse_line_nt_scalar(std::str::from_utf8(line).unwrap()));
            }
        })
    });

    group.bench_function("dispatched", |b| {
        b.iter(|| {
            for line in bytes.split(|&c| c == b'\n') {
                if line.is_empty() { continue; }
                let _ = black_box(M8Record::parse_line_nt(std::str::from_utf8(line).unwrap()));
            }
        })
    });

    group.finish();
}

fn bench_paf_large_batch_streaming(c: &mut Criterion) {
    let batch = generate_paf_batch(500_000);
    let bytes = batch.as_bytes();

    let mut group = c.benchmark_group("paf_large_batch_streaming");
    group.throughput(Throughput::Bytes(bytes.len() as u64));
    group.sample_size(10);

    group.bench_function("scalar", |b| {
        b.iter(|| {
            for line in bytes.split(|&c| c == b'\n') {
                if line.is_empty() { continue; }
                let _ = black_box(PafRecord::parse_line_scalar(std::str::from_utf8(line).unwrap()));
            }
        })
    });

    group.bench_function("dispatched", |b| {
        b.iter(|| {
            for line in bytes.split(|&c| c == b'\n') {
                if line.is_empty() { continue; }
                let _ = black_box(PafRecord::parse_line(std::str::from_utf8(line).unwrap()));
            }
        })
    });

    group.finish();
}


fn generate_fastq_bytes(num_reads: usize, read_len: usize) -> Vec<u8> {
    let mut buffer = Vec::with_capacity(num_reads * (read_len * 2 + 100));
    for i in 0..num_reads {
        let id = format!("read{}", i);
        let seq = "ACGT".repeat(read_len / 4);
        let qual = "IIII".repeat(read_len / 4);

        buffer.extend_from_slice(b"@");
        buffer.extend_from_slice(id.as_bytes());
        buffer.extend_from_slice(b"\n");
        buffer.extend_from_slice(seq.as_bytes());
        buffer.extend_from_slice(b"\n+\n");
        buffer.extend_from_slice(qual.as_bytes());
        buffer.extend_from_slice(b"\n");
    }
    buffer
}

fn bench_read_ingestion_single_end(c: &mut Criterion) {
    let fastq_data = generate_fastq_bytes(200_000, 150);
    let data_len = fastq_data.len() as u64;

    let mut group = c.benchmark_group("read_ingestion_single_end");
    group.throughput(Throughput::Bytes(data_len));
    group.sample_size(10);

    group.bench_function("parse_fastq_to_sequence_record", |b| {
        b.iter(|| {
            let cursor = Cursor::new(black_box(&fastq_data));
            let mut reader = FastqReader::new(cursor);
            let mut count = 0u64;

            while let Some(result) = reader.next() {
                if let Ok(record) = result {
                    let _seq_record: SequenceRecord = record.into(); // includes parse_header
                    count += 1;
                }
            }
            black_box(count)
        })
    });

    group.finish();
}

fn bench_read_ingestion_paired_end_compare_ids(c: &mut Criterion) {
    let mut group = c.benchmark_group("read_ingestion_paired_compare_ids");
    group.sample_size(100);

    let id1 = b"SRR12345678.1 1:N:0:ATCGATCG";
    let id2 = b"SRR12345678.2 2:N:0:ATCGATCG";

    group.bench_function("compare_read_ids_bytes", |b| {
        b.iter(|| {
            // Using the public wrapper for now
            black_box(seqtoid_pipelines::utils::fastx::compare_read_ids(
                black_box(std::str::from_utf8(id1).unwrap()),
                black_box(std::str::from_utf8(id2).unwrap()),
            ))
        })
    });

    group.finish();
}
// ── Register all benchmarks ──────────────────────────────────────────────

criterion_group!(
    benches,
    bench_parse_header,
    bench_paf_line,
    bench_m8_nt_line,
    bench_m8_nr_line,
    bench_m8_nt_large_batch_streaming,
    bench_paf_large_batch_streaming,
    bench_read_ingestion_single_end,
    bench_read_ingestion_paired_end_compare_ids,
);

criterion_main!(benches);