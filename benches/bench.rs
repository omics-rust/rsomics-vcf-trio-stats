use std::hint::black_box;
use std::path::PathBuf;

use criterion::{Criterion, criterion_group, criterion_main};
use rsomics_vcf_trio_stats::{TrioSpec, trio_stats};

fn fixture() -> (PathBuf, PathBuf) {
    let base = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden");
    (base.join("trio.vcf"), base.join("trio.ped"))
}

fn bench_trio_stats(c: &mut Criterion) {
    let (vcf, ped) = fixture();
    c.bench_function("trio_stats_golden", |b| {
        b.iter(|| {
            let t = trio_stats(black_box(&vcf), &TrioSpec::Ped(&ped)).unwrap();
            black_box(t.rows.len())
        });
    });
}

criterion_group!(benches, bench_trio_stats);
criterion_main!(benches);
