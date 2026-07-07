//! Convention normalization demo: fixes QAOA cross-backend fidelity.
//!
//! Demonstrates the CleitonForge finding: quantrs2 uses an Rz sign convention
//! opposite to IBM/Qiskit. Running the same QAOA circuit through native and
//! quantrs2 produces fidelity = 0. After normalization (negating Rz angles),
//! fidelity returns to 1.0.
//!
//! Run:
//!   cargo run --release --example normalize_qaoa -p cforge-cli

use std::f64::consts::PI;

use cforge_backends::{NativeStateVectorBackend, QuantRS2Backend, SimulationBackend};
use cforge_core::{Circuit, GateKind, Operation};
use cforge_parser::{normalize_convention, RzConvention};

fn main() {
    println!("╔══════════════════════════════════════════════════════════════╗");
    println!("║   CleitonForge — Convention Normalization Demo              ║");
    println!("╠══════════════════════════════════════════════════════════════╣");
    println!("║  Finding: quantrs2 Rz convention differs from IBM/Qiskit   ║");
    println!("║  Fix:     normalize_convention(circuit, Reversed, Standard) ║");
    println!("╚══════════════════════════════════════════════════════════════╝");
    println!();

    let circuit = qaoa_maxcut_circuit();

    println!(
        "Circuit: QAOA MaxCut ({} qubits, {} ops)",
        circuit.num_qubits(),
        circuit.operations.len()
    );
    println!("  γ = -3π/4 (cost layer),  β = -π/8 (mixer layer)");
    println!();

    // — Phase 1: raw (un-normalized) ————————————————————————————————————————
    println!("Phase 1 — raw circuit (no normalization)");
    println!("{:-<62}", "");

    let native_sv = run(&NativeStateVectorBackend, &circuit, "native  ");
    let q2_sv_raw = run(&QuantRS2Backend, &circuit, "quantrs2");

    let f_raw = fidelity(&native_sv, &q2_sv_raw);
    println!();
    println!("  Cross-backend fidelity (native ↔ quantrs2): {:.8}", f_raw);
    println!(
        "  Convention divergence: {}",
        if f_raw < 0.01 {
            "CONFIRMED ✗"
        } else {
            "not detected ✓"
        }
    );

    // — Phase 2: normalized ——————————————————————————————————————————————————
    println!();
    println!("Phase 2 — after normalize_convention(Reversed → Standard)");
    println!("{:-<62}", "");

    let normalized = normalize_convention(&circuit, RzConvention::Reversed, RzConvention::Standard);
    println!(
        "  Ops modified: {} Rz-family gates had angles negated",
        count_rz_ops(&circuit)
    );

    let q2_sv_norm = run(&QuantRS2Backend, &normalized, "quantrs2");

    let f_norm = fidelity(&native_sv, &q2_sv_norm);
    println!();
    println!(
        "  Cross-backend fidelity (native ↔ quantrs2): {:.8}",
        f_norm
    );
    println!(
        "  Convention fixed: {}",
        if f_norm > 0.999 {
            "YES ✓"
        } else {
            "NO — unexpected ✗"
        }
    );

    // — Summary ——————————————————————————————————————————————————————————————
    println!();
    println!("Summary");
    println!("{:-<62}", "");
    println!("  Before normalization: fidelity = {:.8}  (broken)", f_raw);
    println!("  After  normalization: fidelity = {:.8}  (fixed) ", f_norm);
    println!();

    if f_norm > 0.999 {
        println!("  Normalization working correctly ✅");
    } else {
        eprintln!("  ERROR: normalization did not restore fidelity");
        std::process::exit(1);
    }
}

// ── Circuit builder ───────────────────────────────────────────────────────────

/// QAOA MaxCut on K₂ (complete graph, 2 nodes) — 1 layer.
///
/// Cost Hamiltonian: H_C = ½(1 − Z₀Z₁) = ½ − ½ Z₀Z₁
/// Mixer: H_B = X₀ + X₁
///
/// Circuit:   H⊗H → CX → Rz(2γ) → CX → Rx(2β)⊗Rx(2β)
fn qaoa_maxcut_circuit() -> Circuit {
    let gamma = -3.0 * PI / 4.0;
    let beta = -PI / 8.0;

    let mut c = Circuit::new(2);
    // Initial superposition
    c.push(Operation::new(GateKind::H, vec![0], vec![]));
    c.push(Operation::new(GateKind::H, vec![1], vec![]));
    // Cost layer: e^{-iγ Z₀Z₁}
    c.push(Operation::new(GateKind::Cx, vec![0, 1], vec![]));
    c.push(Operation::new(GateKind::Rz, vec![1], vec![2.0 * gamma]));
    c.push(Operation::new(GateKind::Cx, vec![0, 1], vec![]));
    // Mixer layer: e^{-iβ (X₀+X₁)}
    c.push(Operation::new(GateKind::Rx, vec![0], vec![2.0 * beta]));
    c.push(Operation::new(GateKind::Rx, vec![1], vec![2.0 * beta]));
    c
}

// ── Helpers ───────────────────────────────────────────────────────────────────

fn run(
    backend: &dyn SimulationBackend,
    circuit: &Circuit,
    label: &str,
) -> Vec<num_complex::Complex64> {
    let result = backend.run(circuit, 0, 0).expect("backend failed");
    let sv = &result.statevector;
    println!(
        "  {label} |ψ⟩ = [{:.4}{:+.4}i, {:.4}{:+.4}i, {:.4}{:+.4}i, {:.4}{:+.4}i]",
        sv[0].re, sv[0].im, sv[1].re, sv[1].im, sv[2].re, sv[2].im, sv[3].re, sv[3].im,
    );
    sv.clone()
}

fn fidelity(a: &[num_complex::Complex64], b: &[num_complex::Complex64]) -> f64 {
    let inner: num_complex::Complex64 = a.iter().zip(b.iter()).map(|(x, y)| x.conj() * y).sum();
    inner.norm_sqr()
}

fn count_rz_ops(circuit: &Circuit) -> usize {
    circuit
        .operations
        .iter()
        .filter(|op| {
            matches!(
                op.kind,
                GateKind::Rz
                    | GateKind::Phase
                    | GateKind::Crz
                    | GateKind::Cp
                    | GateKind::U
                    | GateKind::Cu
            )
        })
        .count()
}
