//! Measurement sampling from a statevector probability distribution.

use std::collections::HashMap;
use num_complex::Complex64;

/// Samples `shots` measurements from the probability distribution |sv|²
/// and returns a map from bitstring (MSB-first, e.g. "01") to count.
///
/// `seed` initialises the deterministic LCG — identical seed + circuit
/// always produces the same counts, enabling reproducible benchmarks.
pub fn sample_counts(sv: &[Complex64], shots: usize, seed: u64) -> HashMap<String, usize> {
    let probs: Vec<f64> = sv.iter().map(|a| a.norm_sqr()).collect();
    let n_qubits = probs.len().trailing_zeros() as usize;

    let mut counts: HashMap<String, usize> = HashMap::new();
    // Use a simple linear-time sampler: build a CDF and binary-search.
    let mut cdf: Vec<f64> = Vec::with_capacity(probs.len());
    let mut acc = 0.0;
    for &p in &probs {
        acc += p;
        cdf.push(acc);
    }

    let mut rng = Lcg64::new(seed);
    for _ in 0..shots {
        let r: f64 = rng.next_f64();
        let idx = cdf.partition_point(|&c| c < r).min(probs.len() - 1);
        let bitstring = format!("{:0>width$b}", idx, width = n_qubits);
        *counts.entry(bitstring).or_insert(0) += 1;
    }
    counts
}

/// Tiny deterministic LCG PRNG — no external randomness dependency.
struct Lcg64 {
    state: u64,
}

impl Lcg64 {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self.state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.state
    }

    fn next_f64(&mut self) -> f64 {
        // Use top 53 bits for double precision.
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}
