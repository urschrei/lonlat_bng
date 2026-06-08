//! Benchmarks comparing the two Transverse Mercator projections.
//!
//! The projection used is selected at compile time by the `karney_tm` feature, so
//! the two are compared across two runs using criterion baselines. The benchmark
//! ids are feature-independent, so the second run reports the delta against the
//! first:
//!
//! ```text
//! cargo bench --bench projection -- --save-baseline redfearn
//! cargo bench --bench projection --features karney_tm -- --baseline redfearn
//! ```
//!
//! `forward_tm`/`inverse_tm` isolate the projection itself (no OSTN15), while the
//! `*_full` variants include the OSTN15 grid shift (forward) and the reverse
//! iteration (inverse), which dominate the end-to-end cost.

use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use rand::SeedableRng;
use rand::distr::{Distribution, Uniform};
use rand::rngs::StdRng;

use lonlat_bng::{convert_etrs89, convert_etrs89_to_ll, convert_osgb36, convert_osgb36_to_ll};

const N: usize = 10_000;

/// A deterministic set of lon/lat points with OSTN15 coverage, so both feature
/// builds benchmark exactly the same inputs.
fn sample_lonlat() -> Vec<(f64, f64)> {
    let mut rng = StdRng::seed_from_u64(0x05D_15u64);
    let lon = Uniform::new(-6.379880, 1.768960).unwrap();
    let lat = Uniform::new(49.871159, 55.811741).unwrap();
    let mut v = Vec::with_capacity(N);
    while v.len() < N {
        let (lo, la) = (lon.sample(&mut rng), lat.sample(&mut rng));
        // Keep only points inside the OSTN15 grid (skip sea/no-coverage points).
        if convert_osgb36(lo, la).is_ok() {
            v.push((lo, la));
        }
    }
    v
}

fn projection_benches(c: &mut Criterion) {
    let ll = sample_lonlat();
    let etrs: Vec<(f64, f64)> = ll
        .iter()
        .map(|&(lo, la)| convert_etrs89(lo, la).unwrap())
        .collect();
    let osgb: Vec<(f64, f64)> = ll
        .iter()
        .map(|&(lo, la)| convert_osgb36(lo, la).unwrap())
        .collect();

    let mut g = c.benchmark_group("projection");

    g.bench_function("forward_tm convert_etrs89", |b| {
        b.iter(|| {
            for &(lo, la) in &ll {
                black_box(convert_etrs89(black_box(lo), black_box(la)).unwrap());
            }
        })
    });
    g.bench_function("inverse_tm convert_etrs89_to_ll", |b| {
        b.iter(|| {
            for &(e, n) in &etrs {
                black_box(convert_etrs89_to_ll(black_box(e), black_box(n)).unwrap());
            }
        })
    });
    g.bench_function("forward_full convert_osgb36", |b| {
        b.iter(|| {
            for &(lo, la) in &ll {
                black_box(convert_osgb36(black_box(lo), black_box(la)).unwrap());
            }
        })
    });
    g.bench_function("inverse_full convert_osgb36_to_ll", |b| {
        b.iter(|| {
            for &(e, n) in &osgb {
                black_box(convert_osgb36_to_ll(black_box(e), black_box(n)).unwrap());
            }
        })
    });

    g.finish();
}

criterion_group!(benches, projection_benches);
criterion_main!(benches);
