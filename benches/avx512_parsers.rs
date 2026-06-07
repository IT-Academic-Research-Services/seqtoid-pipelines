use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use std::fs;

use seqtoid_pipelines::utils::blast::M8Record;
use seqtoid_pipelines::utils::paf::PafRecord;
use seqtoid_pipelines::utils::fastx::parse_header;

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

// ── Register all benchmarks ──────────────────────────────────────────────

criterion_group!(
    benches,
    bench_parse_header,
    bench_paf_line,
    bench_m8_nt_line,
    bench_m8_nr_line,
    bench_m8_nt_large_batch_streaming,
    bench_paf_large_batch_streaming,
);

criterion_main!(benches);