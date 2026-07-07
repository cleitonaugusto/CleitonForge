//! Grover's search algorithm on 3 qubits, comparing both backends.
//!
//! Target state: |101⟩ (state index 5, LSB convention: q0=1, q1=0, q2=1).
//! Optimal iteration count for N=8, M=1: k=2 (sin²(5θ) ≈ 94.8%).
//!
//! Run with:
//!   cargo run --example compare_grover -p cforge-cli

use std::time::Instant;

use cforge_backends::{NativeStateVectorBackend, QuantRS2Backend, SimulationBackend, DEFAULT_SEED};
use cforge_core::{Circuit, GateKind, Operation};
use cforge_metrics::{compute_stats, statevector_fidelity};

/// Build a 3-qubit Grover circuit searching for |101⟩ (state index 5).
///
/// Circuit layout per iteration:
///   - Oracle  : X q1  → CCZ q0,q1,q2  → X q1
///   - Diffuser: H³ → X³ → CCZ → X³ → H³
///
/// CCZ is decomposed as H q2 · CCX q0,q1,q2 · H q2 (phase-kickback).
fn build_grover_circuit() -> Circuit {
    let mut c = Circuit::new(3);

    // Initialize: uniform superposition
    for q in 0..3 {
        c.push(Operation::new(GateKind::H, vec![q], vec![]));
    }

    for _ in 0..2 {
        // ── Oracle: phase-flip |101⟩ (q0=1, q1=0, q2=1) ─────────────────
        // Map |101⟩ → |111⟩, then CCZ, then undo.
        c.push(Operation::new(GateKind::X, vec![1], vec![]));
        c.push(Operation::new(GateKind::H, vec![2], vec![]));
        c.push(Operation::new(GateKind::Ccx, vec![0, 1, 2], vec![]));
        c.push(Operation::new(GateKind::H, vec![2], vec![]));
        c.push(Operation::new(GateKind::X, vec![1], vec![]));

        // ── Diffuser: I − 2|s⟩⟨s|  (inversion about mean) ───────────────
        for q in 0..3 {
            c.push(Operation::new(GateKind::H, vec![q], vec![]));
        }
        for q in 0..3 {
            c.push(Operation::new(GateKind::X, vec![q], vec![]));
        }
        c.push(Operation::new(GateKind::H, vec![2], vec![]));
        c.push(Operation::new(GateKind::Ccx, vec![0, 1, 2], vec![]));
        c.push(Operation::new(GateKind::H, vec![2], vec![]));
        for q in 0..3 {
            c.push(Operation::new(GateKind::X, vec![q], vec![]));
        }
        for q in 0..3 {
            c.push(Operation::new(GateKind::H, vec![q], vec![]));
        }
    }

    c
}

fn main() {
    println!("CleitonForge — Grover search example");
    println!("Target state : |101⟩  (index 5, q0=1 q1=0 q2=1)");
    println!("Qubits       : 3   (N=8 states)");
    println!("Iterations   : 2   (sin²(5θ) ≈ 94.8 % expected)");
    println!();

    let circuit = build_grover_circuit();
    let stats = compute_stats(&circuit);
    println!(
        "Circuit      : {} gates  |  depth {}",
        stats.gate_count, stats.depth
    );
    println!();

    // Reference statevector from the native backend (shots=0 → exact SV).
    let ref_sv = NativeStateVectorBackend
        .run(&circuit, 0, DEFAULT_SEED)
        .expect("native run failed")
        .statevector;

    let backends: [(&str, Box<dyn SimulationBackend>); 2] = [
        ("statevector-native", Box::new(NativeStateVectorBackend)),
        ("statevector-quantrs2", Box::new(QuantRS2Backend)),
    ];

    for (_, backend) in &backends {
        let t0 = Instant::now();
        let result = backend
            .run(&circuit, 1024, DEFAULT_SEED)
            .expect("run failed");
        let elapsed_ms = t0.elapsed().as_secs_f64() * 1000.0;

        // State with highest probability amplitude.
        let (top_idx, top_prob) = result
            .statevector
            .iter()
            .enumerate()
            .map(|(i, a)| (i, a.norm_sqr()))
            .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
            .unwrap();

        // Fidelity relative to native reference.
        let fidelity = statevector_fidelity(&result.statevector, &ref_sv).unwrap_or(1.0);

        // Top 4 bitstrings by shot count.
        let mut counts: Vec<(&String, &usize)> = result.counts.iter().collect();
        counts.sort_by(|a, b| b.1.cmp(a.1));

        println!("Backend  : {}", backend.name());
        println!(
            "  Top state  : |{:03b}⟩  index {}  prob = {:.4} ({:.1} %)",
            top_idx,
            top_idx,
            top_prob,
            top_prob * 100.0
        );
        println!("  Fidelity   : {:.8}", fidelity);
        println!("  Wall time  : {:.3} ms", elapsed_ms);
        println!("  Top counts (1024 shots):");
        for (state, count) in counts.iter().take(4) {
            println!(
                "    |{}⟩  {:>5} shots  ({:.1} %)",
                state,
                count,
                **count as f64 / 1024.0 * 100.0
            );
        }
        println!();
    }

    // Cross-backend fidelity.
    let sv_native = NativeStateVectorBackend
        .run(&circuit, 0, DEFAULT_SEED)
        .expect("native re-run failed")
        .statevector;
    let sv_quantrs2 = QuantRS2Backend
        .run(&circuit, 0, DEFAULT_SEED)
        .expect("quantrs2 re-run failed")
        .statevector;
    let cross = statevector_fidelity(&sv_native, &sv_quantrs2).unwrap_or(0.0);
    println!("Cross-backend fidelity (native vs quantrs2): {:.8}", cross);
    println!(
        "Both backends agree: {}",
        if cross > 0.9999 { "YES ✓" } else { "NO ✗" }
    );
}
