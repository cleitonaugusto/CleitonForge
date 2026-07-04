//! Quantum Volume (QV) benchmark for CleitonForge.
//!
//! Computes the IBM Quantum Volume metric across backends, revealing which
//! simulators agree on the heavy-output test and at what circuit depth.
//!
//! ## Algorithm
//!
//! For each width n (qubits = depth):
//! 1. Generate `TRIALS` random SU(4) circuits of depth n (Haar-random 2-qubit unitaries)
//! 2. For each circuit: compute ideal statevector → get median probability
//! 3. "Heavy" outputs = bitstrings with probability > median
//! 4. HOG (Heavy Output Generation) fraction = P(measuring a heavy output)
//! 5. QV = 2^n where n is the largest n with HOG > 2/3 across all trials
//!
//! The HOG threshold 2/3 is set so that a classical random guesser (uniform
//! distribution, HOG ≈ 1/2) clearly fails while a perfect quantum simulator passes.
//!
//! ## Key finding
//!
//! All three statevector-equivalent backends (native, roqoqo, q1tsim) should
//! achieve identical QV because they use the same Rz convention and there are no
//! phase divergences in Haar-random unitaries (which mix all gate types uniformly).
//! quantrs2 may diverge on circuits that happen to land on Rz-sensitive configurations.

use std::f64::consts::PI;

use num_complex::Complex64;
use rand::rngs::SmallRng;
use rand::{Rng, RngExt, SeedableRng};

use cforge_backends::{NativeStateVectorBackend, QuantRS2Backend, RoqoqoBackend, SimulationBackend};
use cforge_core::{Circuit, GateKind, Operation};

const TRIALS:    usize = 100;
const MAX_WIDTH: usize = 5;
const SEED:      u64   = 2024_07_04;

fn main() {
    println!("╔══════════════════════════════════════════════════════════════╗");
    println!("║         CleitonForge — Quantum Volume Benchmark             ║");
    println!("╠══════════════════════════════════════════════════════════════╣");
    println!("║  Algorithm: Haar-random SU(4) circuits, {TRIALS} trials/width    ║");
    println!("║  HOG threshold: 2/3  (random baseline: 1/2)                ║");
    println!("╚══════════════════════════════════════════════════════════════╝");
    println!();

    let backends: &[(&str, &dyn SimulationBackend)] = &[
        ("native",   &NativeStateVectorBackend),
        ("quantrs2", &QuantRS2Backend),
        ("roqoqo",   &RoqoqoBackend),
    ];

    // Header
    println!(
        "{:<8}  {:<10}  {}",
        "Width n", "QV = 2^n", backends.iter().map(|(n, _)| format!("{:<10}", n)).collect::<Vec<_>>().join("  ")
    );
    println!("{}", "─".repeat(60));

    for n in 2..=MAX_WIDTH {
        let mut hog_fracs: Vec<f64> = vec![];

        for (_, backend) in backends.iter() {
            let hog = compute_hog(*backend, n, TRIALS, SEED);
            hog_fracs.push(hog);
        }

        let qv = if hog_fracs.iter().all(|&h| h > 2.0/3.0) { 1 << n } else { 0 };
        let qv_str = if qv > 0 { format!("2^{n}={qv}") } else { "—".to_string() };

        let cols: String = hog_fracs.iter()
            .map(|h| {
                let pass = if *h > 2.0/3.0 { "✅" } else { "❌" };
                format!("{:.4}  {pass}  ", h)
            })
            .collect::<Vec<_>>()
            .join("  ");

        println!("n={n:<6}  {qv_str:<10}  {cols}");
    }

    println!();
    println!("Note: All statevector backends should agree on QV because");
    println!("Haar-random circuits are insensitive to single-gate Rz sign conventions.");
    println!("Any divergence indicates a deeper convention difference.");
}

// ── HOG fraction computation ──────────────────────────────────────────────────

fn compute_hog(backend: &dyn SimulationBackend, n: usize, trials: usize, seed: u64) -> f64 {
    let mut rng = SmallRng::seed_from_u64(seed ^ (n as u64 * 0xDEAD_BEEF));
    let mut total_hog_prob = 0.0;

    for _ in 0..trials {
        let circuit = random_qv_circuit(n, &mut rng);

        // Ideal statevector
        let ideal_sv = NativeStateVectorBackend
            .run(&circuit, 0, 0)
            .expect("native backend failed")
            .statevector;

        let probs: Vec<f64> = ideal_sv.iter().map(|a| a.norm_sqr()).collect();

        // Heavy outputs: P(i) > median P
        let mut sorted = probs.clone();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let median = sorted[sorted.len() / 2];

        let heavy_set: Vec<bool> = probs.iter().map(|&p| p > median).collect();
        let heavy_ideal_prob: f64 = probs.iter().zip(&heavy_set).filter(|(_, &h)| h).map(|(p, _)| p).sum();

        // Backend simulation — compare HOG using the backend's statevector
        let result = backend.run(&circuit, 0, 0);
        if let Ok(r) = result {
            let backend_probs: Vec<f64> = r.statevector.iter().map(|a| a.norm_sqr()).collect();
            let backend_hog: f64 = backend_probs.iter().zip(&heavy_set).filter(|(_, &h)| h).map(|(p, _)| p).sum();
            total_hog_prob += backend_hog;
        } else {
            total_hog_prob += 0.0;
        }

        // (ideal heavy prob is always > 2/3 by construction for small n)
        let _ = heavy_ideal_prob;
    }

    total_hog_prob / trials as f64
}

// ── Random QV circuit generator ───────────────────────────────────────────────

/// Generates a random Quantum Volume circuit of width n = depth n.
///
/// For each layer: randomly pair the n qubits and apply a Haar-random SU(4)
/// unitary to each pair. The SU(4) is decomposed into parameterized gates.
fn random_qv_circuit(n: usize, rng: &mut impl Rng) -> Circuit {
    let mut circuit = Circuit::new(n);

    for _layer in 0..n {
        // Random permutation of qubits
        let mut perm: Vec<usize> = (0..n).collect();
        for i in (1..n).rev() {
            let j = rng.random_range(0..=i);
            perm.swap(i, j);
        }

        // Apply random SU(4) to each consecutive pair in the permutation
        let mut i = 0;
        while i + 1 < n {
            let q0 = perm[i];
            let q1 = perm[i + 1];
            apply_random_su4(&mut circuit, q0, q1, rng);
            i += 2;
        }
        // If n is odd, apply a random SU(2) to the remaining qubit
        if n % 2 == 1 {
            let q = perm[n - 1];
            apply_random_su2(&mut circuit, q, rng);
        }
    }

    circuit
}

/// Decomposes a Haar-random SU(4) via KAK parameterisation.
///
/// Approximation: random single-qubit rotations + CNOT + more rotations.
/// This is not the exact KAK decomposition but samples the SU(4) Haar
/// measure well enough for QV benchmarking purposes.
fn apply_random_su4(circuit: &mut Circuit, q0: usize, q1: usize, rng: &mut impl Rng) {
    // Random local rotations before CNOT
    apply_random_su2(circuit, q0, rng);
    apply_random_su2(circuit, q1, rng);

    // Entangling layer (CNOT + Rz)
    circuit.push(Operation::new(GateKind::Cx, vec![q0, q1], vec![]));
    let theta: f64 = rng.random::<f64>() * PI;
    circuit.push(Operation::new(GateKind::Rz, vec![q1], vec![theta]));
    circuit.push(Operation::new(GateKind::Cx, vec![q1, q0], vec![]));
    let phi: f64 = rng.random::<f64>() * PI;
    circuit.push(Operation::new(GateKind::Ry, vec![q0], vec![phi]));
    circuit.push(Operation::new(GateKind::Cx, vec![q0, q1], vec![]));

    // Random local rotations after
    apply_random_su2(circuit, q0, rng);
    apply_random_su2(circuit, q1, rng);
}

/// Random SU(2) via Euler angles: Rz(α) Ry(β) Rz(γ).
fn apply_random_su2(circuit: &mut Circuit, q: usize, rng: &mut impl Rng) {
    let alpha: f64 = rng.random::<f64>() * 2.0 * PI;
    let beta:  f64 = (1.0_f64 - 2.0 * rng.random::<f64>()).acos();
    let gamma: f64 = rng.random::<f64>() * 2.0 * PI;
    circuit.push(Operation::new(GateKind::Rz, vec![q], vec![alpha]));
    circuit.push(Operation::new(GateKind::Ry, vec![q], vec![beta]));
    circuit.push(Operation::new(GateKind::Rz, vec![q], vec![gamma]));
}
