//! Quantum Volume (QV) benchmark for CleitonForge.
//!
//! Computes the IBM Quantum Volume metric across backends, revealing which
//! simulators agree on the heavy-output test and at what circuit depth.
//!
//! ## Algorithm
//!
//! For each width n (qubits = depth):
//! 1. Generate `TRIALS` random SU(4) circuits of depth n via KAK decomposition
//! 2. For each circuit: compute ideal statevector → get median probability
//! 3. "Heavy" outputs = bitstrings with probability > median
//! 4. HOG (Heavy Output Generation) fraction = P(measuring a heavy output)
//! 5. QV = 2^n where n is the largest n with HOG > 2/3 across all trials
//!
//! The HOG threshold 2/3 is set so that a classical random guesser (uniform
//! distribution, HOG ≈ 1/2) clearly fails while a perfect quantum simulator passes.
//!
//! ## SU(4) circuit generation
//!
//! Each 2-qubit block uses the KAK decomposition:
//!   U = (A₁⊗A₂) · exp(-i(c₁·XX + c₂·YY + c₃·ZZ)) · (B₁⊗B₂)
//! The SU(2) factors are Haar-exact (Euler angles with correct spherical measure).
//! The Weyl interaction uses three parameters (c₁,c₂,c₃) sampled via order statistics
//! of three iid U[0,π/4] draws, giving the uniform distribution over the ordered
//! simplex {0≤c₃≤c₂≤c₁≤π/4} — covers all of SU(4) but with non-Haar weight.
//!
//! ## Key finding
//!
//! QV heavy-output scores are structurally blind to Rz sign bugs when circuits
//! are built from Rz + real-Ry + real-CX: swapping Rz(θ)→Rz(-θ) complex-conjugates
//! the entire unitary (U_wrong = U*_ref), leaving all computational-basis probabilities
//! exactly unchanged. This is an exact blindness, not a statistical average.

use std::f64::consts::PI;

use rand::rngs::SmallRng;
use rand::{Rng, RngExt, SeedableRng};

use cforge_backends::{
    NativeStateVectorBackend, Q1tSimBackend, QuantRS2Backend, RoqoqoBackend, SimulationBackend,
};
use cforge_core::{Circuit, GateKind, Operation};

const TRIALS: usize = 100;
const MAX_WIDTH: usize = 5;
const SEED: u64 = 2024_07_04;

fn main() {
    println!("╔══════════════════════════════════════════════════════════════╗");
    println!("║         CleitonForge — Quantum Volume Benchmark             ║");
    println!("╠══════════════════════════════════════════════════════════════╣");
    println!("║  Algorithm: KAK SU(4) circuits (uniform Weyl), {TRIALS} trials/width ║");
    println!("║  HOG threshold: 2/3  (random baseline: 1/2)                ║");
    println!("╚══════════════════════════════════════════════════════════════╝");
    println!();

    let backends: &[(&str, &dyn SimulationBackend)] = &[
        ("native", &NativeStateVectorBackend),
        ("quantrs2", &QuantRS2Backend),
        ("roqoqo", &RoqoqoBackend),
        ("q1tsim", &Q1tSimBackend),
    ];

    // Header
    println!(
        "{:<8}  {:<10}  {}",
        "Width n",
        "QV = 2^n",
        backends
            .iter()
            .map(|(n, _)| format!("{:<10}", n))
            .collect::<Vec<_>>()
            .join("  ")
    );
    println!("{}", "─".repeat(60));

    for n in 2..=MAX_WIDTH {
        let mut hog_fracs: Vec<f64> = vec![];

        for (_, backend) in backends.iter() {
            let hog = compute_hog(*backend, n, TRIALS, SEED);
            hog_fracs.push(hog);
        }

        let qv = if hog_fracs.iter().all(|&h| h > 2.0 / 3.0) {
            1 << n
        } else {
            0
        };
        let qv_str = if qv > 0 {
            format!("2^{n}={qv}")
        } else {
            "—".to_string()
        };

        let cols: String = hog_fracs
            .iter()
            .map(|h| {
                let pass = if *h > 2.0 / 3.0 { "✅" } else { "❌" };
                format!("{:.4}  {pass}  ", h)
            })
            .collect::<Vec<_>>()
            .join("  ");

        println!("n={n:<6}  {qv_str:<10}  {cols}");
    }

    println!();
    println!("Note: QV heavy-output scores are blind to Rz sign bugs (U_wrong = U*_ref");
    println!("gives identical probabilities). Use phase-sensitive conformance tests");
    println!("(e.g. Rz(π/2)|+⟩ ~ (|0⟩+i|1⟩)/√2) for cross-backend sign verification.");
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

        // Heavy outputs: P(i) > median P (lower median for even-length arrays,
        // matching IBM's QV definition: median = average of two central elements).
        let mut sorted = probs.clone();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let mid = sorted.len() / 2;
        let median = (sorted[mid - 1] + sorted[mid]) / 2.0;

        let heavy_set: Vec<bool> = probs.iter().map(|&p| p > median).collect();

        // Backend simulation — compare HOG using the backend's statevector
        let result = backend.run(&circuit, 0, 0);
        if let Ok(r) = result {
            assert_eq!(
                r.statevector.len(),
                probs.len(),
                "backend returned statevector of length {} but expected {}",
                r.statevector.len(),
                probs.len()
            );
            let backend_hog: f64 = r
                .statevector
                .iter()
                .zip(&heavy_set)
                .filter(|(_, &h)| h)
                .map(|(a, _)| a.norm_sqr())
                .sum();
            total_hog_prob += backend_hog;
        }
    }

    total_hog_prob / trials as f64
}

// ── Random QV circuit generator ───────────────────────────────────────────────

/// Generates a random Quantum Volume circuit of width n = depth n.
///
/// For each layer: randomly pair the n qubits and apply a random SU(4) block
/// via the KAK decomposition with uniform Weyl chamber sampling.
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

/// Applies a random SU(4) via the KAK (Cartan) decomposition:
///   U = (A₁⊗A₂) · exp(-i(c₁·XX + c₂·YY + c₃·ZZ)) · (B₁⊗B₂)
///
/// The SU(2) factors (A₁,A₂,B₁,B₂) are Haar-exact via Euler angles with the
/// correct spherical measure. The Weyl parameters (c₁,c₂,c₃) are the order
/// statistics of three iid U[0,π/4] draws — uniform over {0≤c₃≤c₂≤c₁≤π/4},
/// covering all 15 real dimensions of SU(4) with a non-Haar interaction weight.
///
/// Since XX, YY, ZZ mutually commute, the interaction factors as:
///   exp(-i·c₁·XX) · exp(-i·c₂·YY) · exp(-i·c₃·ZZ)
/// Each ZZ-type block uses CX·Rz(2c)·CX = exp(-i·c·ZZ) (from Rz(θ)=diag(e^{-iθ/2},e^{+iθ/2})).
/// The negative sign is consistent across all three blocks; QV scores are unaffected.
fn apply_random_su4(circuit: &mut Circuit, q0: usize, q1: usize, rng: &mut impl Rng) {
    apply_random_su2(circuit, q0, rng);
    apply_random_su2(circuit, q1, rng);

    // Weyl chamber: c₁ ≥ c₂ ≥ c₃ ≥ 0, c₁ ∈ [0, π/4].
    // Use order statistics of 3 iid U[0,π/4] draws — gives the uniform
    // distribution over the ordered simplex {0≤c₃≤c₂≤c₁≤π/4}.
    let mut raw = [
        rng.random::<f64>() * PI / 4.0,
        rng.random::<f64>() * PI / 4.0,
        rng.random::<f64>() * PI / 4.0,
    ];
    raw.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
    let [c1, c2, c3] = raw;

    // exp(-i·c₃·ZZ): CX · Rz(2c₃) · CX
    // Note: with Rz(θ)=diag(e^{-iθ/2},e^{+iθ/2}), this gives diag(e^{-ic},e^{+ic},e^{+ic},e^{-ic})
    // = exp(-i·c·ZZ). The negative sign is consistent across all three blocks; QV scores
    // (which depend only on |amplitude|²) are unaffected.
    circuit.push(Operation::new(GateKind::Cx, vec![q0, q1], vec![]));
    circuit.push(Operation::new(GateKind::Rz, vec![q1], vec![2.0 * c3]));
    circuit.push(Operation::new(GateKind::Cx, vec![q0, q1], vec![]));

    // exp(-i·c₂·YY): Rx(π/2)·Z·Rx(-π/2) = -Y, so
    //   (Rx(π/2)⊗Rx(π/2)) · exp(-i·c₂·ZZ) · (Rx(-π/2)⊗Rx(-π/2)) = exp(-i·c₂·YY)
    let hpi = PI / 2.0;
    circuit.push(Operation::new(GateKind::Rx, vec![q0], vec![hpi]));
    circuit.push(Operation::new(GateKind::Rx, vec![q1], vec![hpi]));
    circuit.push(Operation::new(GateKind::Cx, vec![q0, q1], vec![]));
    circuit.push(Operation::new(GateKind::Rz, vec![q1], vec![2.0 * c2]));
    circuit.push(Operation::new(GateKind::Cx, vec![q0, q1], vec![]));
    circuit.push(Operation::new(GateKind::Rx, vec![q0], vec![-hpi]));
    circuit.push(Operation::new(GateKind::Rx, vec![q1], vec![-hpi]));

    // exp(-i·c₁·XX): H·Z·H = X, so
    //   (H⊗H) · exp(-i·c₁·ZZ) · (H⊗H) = exp(-i·c₁·XX)
    circuit.push(Operation::new(GateKind::H, vec![q0], vec![]));
    circuit.push(Operation::new(GateKind::H, vec![q1], vec![]));
    circuit.push(Operation::new(GateKind::Cx, vec![q0, q1], vec![]));
    circuit.push(Operation::new(GateKind::Rz, vec![q1], vec![2.0 * c1]));
    circuit.push(Operation::new(GateKind::Cx, vec![q0, q1], vec![]));
    circuit.push(Operation::new(GateKind::H, vec![q0], vec![]));
    circuit.push(Operation::new(GateKind::H, vec![q1], vec![]));

    apply_random_su2(circuit, q0, rng);
    apply_random_su2(circuit, q1, rng);
}

/// Random SU(2) via Euler angles: Rz(α) Ry(β) Rz(γ).
fn apply_random_su2(circuit: &mut Circuit, q: usize, rng: &mut impl Rng) {
    let alpha: f64 = rng.random::<f64>() * 2.0 * PI;
    let beta: f64 = (1.0_f64 - 2.0 * rng.random::<f64>()).acos();
    let gamma: f64 = rng.random::<f64>() * 2.0 * PI;
    circuit.push(Operation::new(GateKind::Rz, vec![q], vec![alpha]));
    circuit.push(Operation::new(GateKind::Ry, vec![q], vec![beta]));
    circuit.push(Operation::new(GateKind::Rz, vec![q], vec![gamma]));
}
