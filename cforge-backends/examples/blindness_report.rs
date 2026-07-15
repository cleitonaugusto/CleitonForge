//! Empirical benchmark-blindness report.
//!
//! Runs every standard benchmark family against a correct simulator and
//! against two deliberately buggy ones — a full conjugation bug (every
//! gate implemented as U*) and the quantrs2-v0.2.0-class rotation-sign
//! bug — and prints the scores side by side.
//!
//! Expected result, per the conjugation-invariance theorem: QV-HOG,
//! XEB, mirror circuits and even full Quantum Phase Estimation return
//! identical scores for the correct and the bugged simulators. Only the
//! amplitude-level QGCS checks (and, for the partial bug, a crafted
//! 4-gate witness) separate them.
//!
//! Run with: cargo run -p cforge-backends --example blindness_report

use std::f64::consts::PI;

use cforge_backends::{
    certify, heavy_output_probability, inverse_circuit, mirror_survival, random_qv_circuit,
    xeb_score_exact, xeb_score_sampled, ConjugatedStateVectorBackend, ConjugationScope,
    NativeStateVectorBackend, SimulationBackend,
};
use cforge_core::{Circuit, GateKind, Operation};

const SHOTS: usize = 100_000;
const N_CIRCUITS: u64 = 20;

fn main() {
    let ideal = NativeStateVectorBackend;
    let conj_all = ConjugatedStateVectorBackend::new(ConjugationScope::AllGates);
    let conj_rz = ConjugatedStateVectorBackend::new(ConjugationScope::RotationSignOnly);

    println!("════════════════════════════════════════════════════════════════");
    println!(" Benchmark blindness report — CleitonForge");
    println!(" Correct simulator vs U*-conjugated vs Rz-sign bug (quantrs2 class)");
    println!("════════════════════════════════════════════════════════════════\n");

    // ── 1. Quantum Volume HOG ────────────────────────────────────────
    println!("1. Quantum Volume — heavy-output-generation probability");
    println!("   ({} random 5q×5 Rz/Ry/CX circuits, {} shots each)\n", N_CIRCUITS, SHOTS);
    println!("   {:<8} {:>12} {:>12} {:>12}", "seed", "correct", "U* all", "Rz-sign");
    let mut max_hog_gap: f64 = 0.0;
    for seed in 0..N_CIRCUITS {
        let c = random_qv_circuit(5, 5, seed);
        let h0 = heavy_output_probability(&c, &ideal, &ideal, SHOTS, seed).unwrap();
        let h1 = heavy_output_probability(&c, &ideal, &conj_all, SHOTS, seed).unwrap();
        let h2 = heavy_output_probability(&c, &ideal, &conj_rz, SHOTS, seed).unwrap();
        max_hog_gap = max_hog_gap.max((h0 - h1).abs()).max((h0 - h2).abs());
        if seed < 5 {
            println!("   {:<8} {:>12.6} {:>12.6} {:>12.6}", seed, h0, h1, h2);
        }
    }
    println!("   …");
    println!("   max |Δ HOG| across all {} circuits: {:.2e}\n", N_CIRCUITS, max_hog_gap);

    // ── 2. Linear XEB ────────────────────────────────────────────────
    println!("2. Linear cross-entropy benchmarking (Google supremacy metric)\n");
    println!(
        "   {:<8} {:>14} {:>14} {:>14}",
        "seed", "correct", "U* all", "Rz-sign"
    );
    let mut max_xeb_gap: f64 = 0.0;
    for seed in 0..N_CIRCUITS {
        let c = random_qv_circuit(5, 5, seed);
        let f0 = xeb_score_exact(&c, &ideal, &ideal).unwrap();
        let f1 = xeb_score_exact(&c, &ideal, &conj_all).unwrap();
        let f2 = xeb_score_exact(&c, &ideal, &conj_rz).unwrap();
        max_xeb_gap = max_xeb_gap.max((f0 - f1).abs()).max((f0 - f2).abs());
        if seed < 5 {
            println!("   {:<8} {:>14.9} {:>14.9} {:>14.9}", seed, f0, f1, f2);
        }
    }
    println!("   …");
    println!("   max |Δ F_XEB| (exact) across all circuits: {:.2e}", max_xeb_gap);
    let c = random_qv_circuit(5, 5, 3);
    let fs0 = xeb_score_sampled(&c, &ideal, &ideal, SHOTS, 99).unwrap();
    let fs1 = xeb_score_sampled(&c, &ideal, &conj_all, SHOTS, 99).unwrap();
    println!(
        "   sampled ({} shots): correct = {:.6}, U* = {:.6}, Δ = {:.2e}\n",
        SHOTS,
        fs0,
        fs1,
        (fs0 - fs1).abs()
    );

    // ── 3. Mirror circuits ───────────────────────────────────────────
    println!("3. Mirror circuits (Sandia) — survival probability P(|0…0⟩)\n");
    let c = random_qv_circuit(5, 5, 7);
    let m0 = mirror_survival(&c, &ideal).unwrap();
    let m1 = mirror_survival(&c, &conj_all).unwrap();
    let m2 = mirror_survival(&c, &conj_rz).unwrap();
    println!("   correct: {m0:.12}");
    println!("   U* all:  {m1:.12}");
    println!("   Rz-sign: {m2:.12}\n");

    // ── 4. Quantum Phase Estimation ──────────────────────────────────
    println!("4. Quantum Phase Estimation — the strongest case");
    println!("   Estimating φ = 5/8 of Phase(2πφ), 3 counting qubits.");
    println!("   Correct answer: counting register = |101⟩ = 5.\n");
    let qpe = qpe_circuit(3, 5.0 / 8.0);
    for backend in [
        &ideal as &dyn SimulationBackend,
        &conj_all as &dyn SimulationBackend,
        &conj_rz as &dyn SimulationBackend,
    ] {
        let r = backend.run(&qpe, 0, 0).unwrap();
        let probs = r.probabilities();
        let (best, p) = probs
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
            .unwrap();
        let counting = best & 0b111; // low 3 bits = counting register
        println!(
            "   {:<26} peak index {:>2} → counting = {} (P = {:.9})",
            backend.name(),
            best,
            counting,
            p
        );
    }
    println!("\n   All three report the CORRECT phase. A simulator in which every");
    println!("   gate is wrong executes textbook QPE flawlessly — global");
    println!("   conjugation is a gauge symmetry of measurement statistics.\n");

    // ── 5. What DOES see the bugs ────────────────────────────────────
    println!("5. Detection — what the benchmarks miss\n");

    // 5a. Crafted witness for the partial bug.
    let theta = 1.0;
    let mut w = Circuit::new(1);
    w.push(Operation::new(GateKind::H, vec![0], vec![]));
    w.push(Operation::new(GateKind::Rz, vec![0], vec![theta]));
    w.push(Operation::new(GateKind::S, vec![0], vec![]));
    w.push(Operation::new(GateKind::H, vec![0], vec![]));
    let p_good = ideal.run(&w, 0, 0).unwrap().probabilities()[0];
    let p_bad = conj_rz.run(&w, 0, 0).unwrap().probabilities()[0];
    println!("   4-gate witness H·Rz(1)·S·H for the Rz-sign bug:");
    println!(
        "     P(0) correct = {:.6}, buggy = {:.6}, Δp = {:.6}  ← visible in counts",
        p_good,
        p_bad,
        (p_good - p_bad).abs()
    );
    let p_bad_all = conj_all.run(&w, 0, 0).unwrap().probabilities()[0];
    println!(
        "     (same witness vs U*-all bug: Δp = {:.2e} — full conjugation stays invisible)\n",
        (p_good - p_bad_all).abs()
    );

    // 5b. QGCS amplitude-level checks.
    for (label, backend) in [
        ("U*-conjugated", &conj_all as &dyn SimulationBackend),
        ("Rz-sign bug", &conj_rz as &dyn SimulationBackend),
    ] {
        let results = certify(backend);
        let failed: Vec<&str> = results
            .iter()
            .filter(|r| !r.passed() && !r.skipped())
            .map(|r| r.name)
            .collect();
        let total = results.len();
        println!(
            "   QGCS certify({label}): {}/{} checks FAIL → {}",
            failed.len(),
            total,
            failed.join(", ")
        );
    }

    // Sanity: the correct backend passes everything.
    let clean = certify(&ideal);
    let clean_fails = clean.iter().filter(|r| !r.passed() && !r.skipped()).count();
    println!("   QGCS certify(correct): {clean_fails} failures\n");

    println!("════════════════════════════════════════════════════════════════");
    println!(" Conclusion: QV, XEB, mirror and QPE cannot distinguish a correct");
    println!(" simulator from a conjugated one — provably. Amplitude-level");
    println!(" conformance (QGCS) is the only test here that can.");
    println!("════════════════════════════════════════════════════════════════");
}

// ── QPE construction (little-endian: qubit 0 = LSB) ──────────────────

/// QFT over `m` qubits, |x⟩ → 2^{−m/2} Σ_y e^{2πi·xy/2^m}|y⟩.
fn qft_circuit(m: usize) -> Circuit {
    let mut c = Circuit::new(m);
    for j in (0..m).rev() {
        c.push(Operation::new(GateKind::H, vec![j], vec![]));
        for k in (0..j).rev() {
            let angle = PI / (1 << (j - k)) as f64;
            c.push(Operation::new(GateKind::Cp, vec![k, j], vec![angle]));
        }
    }
    for j in 0..m / 2 {
        c.push(Operation::new(GateKind::Swap, vec![j, m - 1 - j], vec![]));
    }
    c
}

/// Textbook QPE for U = Phase(2πφ) with `m` counting qubits (0..m) and
/// the eigenstate |1⟩ prepared on qubit `m`.
fn qpe_circuit(m: usize, phi: f64) -> Circuit {
    let mut c = Circuit::new(m + 1);
    c.push(Operation::new(GateKind::X, vec![m], vec![]));
    for j in 0..m {
        c.push(Operation::new(GateKind::H, vec![j], vec![]));
    }
    for j in 0..m {
        let angle = 2.0 * PI * phi * (1 << j) as f64;
        c.push(Operation::new(GateKind::Cp, vec![j, m], vec![angle]));
    }
    // Inverse QFT on the counting register, obtained by exact IR-level
    // inversion of the forward QFT.
    let qft_inv = inverse_circuit(&qft_circuit(m));
    for op in qft_inv.operations {
        c.push(op);
    }
    c
}
