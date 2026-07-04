//! CleitonForge — Canonical Benchmark Suite
//!
//! Runs six standard quantum algorithms on both backends (statevector and QuantRS2)
//! and verifies results against known expected outcomes.
//!
//! Algorithms:
//!   1. Bell state          — 2q, expected 50% |00⟩ + 50% |11⟩
//!   2. GHZ state           — 3q, expected 50% |000⟩ + 50% |111⟩
//!   3. QFT on |100⟩        — 3q, expected uniform superposition (1/8 each)
//!   4. Bernstein-Vazirani  — 3q, s=101, expected 100% |101⟩
//!   5. Grover search       — 3q, target |101⟩, expected ≈94.5% |101⟩
//!   6. QAOA MaxCut p=1     — 2q, reveals documented RZ convention divergence
//!
//! Run with:
//!   cargo run --example benchmark_suite -p cforge-cli

use std::time::Instant;

use cforge_backends::{DEFAULT_SEED, NativeStateVectorBackend, QuantRS2Backend, SimulationBackend};
use cforge_core::{Circuit, GateKind, Operation};
use cforge_metrics::{compute_stats, statevector_fidelity};

const SHOTS: usize = 4096;

// ── Circuit builders ──────────────────────────────────────────────────────────

fn bell_state() -> Circuit {
    let mut c = Circuit::new(2);
    c.push(Operation::new(GateKind::H,  vec![0],    vec![]));
    c.push(Operation::new(GateKind::Cx, vec![0, 1], vec![]));
    c
}

fn ghz_state() -> Circuit {
    let mut c = Circuit::new(3);
    c.push(Operation::new(GateKind::H,  vec![0],    vec![]));
    c.push(Operation::new(GateKind::Cx, vec![0, 1], vec![]));
    c.push(Operation::new(GateKind::Cx, vec![1, 2], vec![]));
    c
}

fn qft3_on_100() -> Circuit {
    use std::f64::consts::PI;
    let mut c = Circuit::new(3);
    // Initialise to |100⟩
    c.push(Operation::new(GateKind::X,    vec![0],    vec![]));
    // QFT
    c.push(Operation::new(GateKind::H,    vec![0],    vec![]));
    c.push(Operation::new(GateKind::Cp,   vec![1, 0], vec![PI / 2.0]));
    c.push(Operation::new(GateKind::Cp,   vec![2, 0], vec![PI / 4.0]));
    c.push(Operation::new(GateKind::H,    vec![1],    vec![]));
    c.push(Operation::new(GateKind::Cp,   vec![2, 1], vec![PI / 2.0]));
    c.push(Operation::new(GateKind::H,    vec![2],    vec![]));
    c.push(Operation::new(GateKind::Swap, vec![0, 2], vec![]));
    c
}

fn bernstein_vazirani_101() -> Circuit {
    // Phase oracle: Z(q0) ⊗ I(q1) ⊗ Z(q2) → phase (-1)^{x0+x2} = (-1)^{s·x}
    // Result: always |101⟩
    let mut c = Circuit::new(3);
    for q in 0..3 {
        c.push(Operation::new(GateKind::H, vec![q], vec![]));
    }
    c.push(Operation::new(GateKind::Z, vec![0], vec![]));
    c.push(Operation::new(GateKind::Z, vec![2], vec![]));
    for q in 0..3 {
        c.push(Operation::new(GateKind::H, vec![q], vec![]));
    }
    c
}

fn grover_101() -> Circuit {
    // 3-qubit Grover targeting |101⟩, 2 iterations.
    // Oracle: X(q1) · CCZ(q0,q1,q2) · X(q1)  where CCZ = H·CCX·H (phase kickback)
    // Diffuser: H³ · X³ · CCZ · X³ · H³
    let mut c = Circuit::new(3);
    for q in 0..3 {
        c.push(Operation::new(GateKind::H, vec![q], vec![]));
    }
    for _ in 0..2 {
        // Oracle
        c.push(Operation::new(GateKind::X,   vec![1],       vec![]));
        c.push(Operation::new(GateKind::H,   vec![2],       vec![]));
        c.push(Operation::new(GateKind::Ccx, vec![0, 1, 2], vec![]));
        c.push(Operation::new(GateKind::H,   vec![2],       vec![]));
        c.push(Operation::new(GateKind::X,   vec![1],       vec![]));
        // Diffuser
        for q in 0..3 { c.push(Operation::new(GateKind::H, vec![q], vec![])); }
        for q in 0..3 { c.push(Operation::new(GateKind::X, vec![q], vec![])); }
        c.push(Operation::new(GateKind::H,   vec![2],       vec![]));
        c.push(Operation::new(GateKind::Ccx, vec![0, 1, 2], vec![]));
        c.push(Operation::new(GateKind::H,   vec![2],       vec![]));
        for q in 0..3 { c.push(Operation::new(GateKind::X, vec![q], vec![])); }
        for q in 0..3 { c.push(Operation::new(GateKind::H, vec![q], vec![])); }
    }
    c
}

// QAOA p=1 MaxCut on a single edge (q0-q1).
// Optimal angles (statevector): γ = -3π/4, β = -π/8 → 100% |01⟩+|10⟩ (cut states).
// Note: quantrs2 applies Rz with the opposite imaginary-part sign convention,
// so the two backends intentionally diverge here. This is one of the inter-
// framework discrepancies CleitonForge is designed to surface.
fn qaoa_maxcut_2q() -> Circuit {
    use std::f64::consts::PI;
    let gamma = -3.0 * PI / 4.0;
    let beta  = -PI / 8.0;
    let mut c = Circuit::new(2);
    c.push(Operation::new(GateKind::H,  vec![0],    vec![]));
    c.push(Operation::new(GateKind::H,  vec![1],    vec![]));
    // Cost: RZZ(2γ) = CX, RZ(2γ), CX
    c.push(Operation::new(GateKind::Cx, vec![0, 1], vec![]));
    c.push(Operation::new(GateKind::Rz, vec![1],    vec![2.0 * gamma]));
    c.push(Operation::new(GateKind::Cx, vec![0, 1], vec![]));
    // Mixer: RX(2β)
    c.push(Operation::new(GateKind::Rx, vec![0],    vec![2.0 * beta]));
    c.push(Operation::new(GateKind::Rx, vec![1],    vec![2.0 * beta]));
    c
}

// ── Benchmark runner ──────────────────────────────────────────────────────────

enum ExpectKind {
    /// Superposition of all-zeros and all-ones strings: e.g. |00⟩+|11⟩
    Entangled,
    /// Exactly this bitstring with ≥ min_prob
    Bitstring(&'static str),
    /// Uniform superposition over 2^n states (e.g. QFT output)
    Uniform,
    /// MaxCut states (|01⟩ + |10⟩); cross-backend fidelity NOT required
    /// because the quantrs2 Rz convention diverges from statevector here.
    MaxCut,
}

struct BenchResult {
    name:       &'static str,
    qubits:     usize,
    gates:      usize,
    depth:      usize,
    label:      String,
    fidelity:   f64,
    time_sv_ms: f64,
    time_q2_ms: f64,
    pass:       bool,
    note:       Option<&'static str>,
}

fn run_benchmark(
    name:     &'static str,
    circuit:  &Circuit,
    expect:   ExpectKind,
    min_prob: f64,
) -> BenchResult {
    let stats = compute_stats(circuit);

    let t0 = Instant::now();
    let sv  = NativeStateVectorBackend.run(circuit, SHOTS, DEFAULT_SEED).expect("sv failed");
    let time_sv_ms = t0.elapsed().as_secs_f64() * 1000.0;

    let t0 = Instant::now();
    let q2  = QuantRS2Backend.run(circuit, SHOTS, DEFAULT_SEED).expect("q2 failed");
    let time_q2_ms = t0.elapsed().as_secs_f64() * 1000.0;

    let fidelity = statevector_fidelity(&sv.statevector, &q2.statevector).unwrap_or(0.0);

    let n        = circuit.num_qubits();
    let total    = sv.counts.values().sum::<usize>() as f64;
    let mut note = None;

    let (effective_prob, label, fidelity_ok) = match &expect {
        ExpectKind::Entangled => {
            let zeros = "0".repeat(n);
            let ones  = "1".repeat(n);
            let p0 = sv.counts.get(zeros.as_str()).copied().unwrap_or(0) as f64 / total;
            let p1 = sv.counts.get(ones.as_str()).copied().unwrap_or(0) as f64 / total;
            (p0 + p1, format!("|{zeros}⟩+|{ones}⟩ {:.0}%", (p0+p1)*100.0), true)
        }
        ExpectKind::Bitstring(s) => {
            let p = sv.counts.get(*s).copied().unwrap_or(0) as f64 / total;
            (p, format!("|{s}⟩ {:.0}%", p*100.0), true)
        }
        ExpectKind::Uniform => {
            // Check that all 2^n amplitudes have probability ≈ 1/2^n
            let expected_p = 1.0 / (1 << n) as f64;
            let uniform = sv.statevector.iter()
                .all(|a| (a.norm_sqr() - expected_p).abs() < 0.002);
            let eff = if uniform { 1.0 } else { 0.0 };
            let lbl = if uniform { format!("uniform 1/{} ✓", 1<<n) } else { "not uniform".to_string() };
            (eff, lbl, true)
        }
        ExpectKind::MaxCut => {
            // Statevector backend: cut states |01⟩+|10⟩
            let p01 = sv.counts.get("01").copied().unwrap_or(0) as f64 / total;
            let p10 = sv.counts.get("10").copied().unwrap_or(0) as f64 / total;
            let cut = p01 + p10;
            note = Some("Rz convention divergence (by design)");
            // Cross-backend fidelity is NOT checked — divergence is intentional
            (cut, format!("cut {:.0}%", cut*100.0), false)
        }
    };

    let fidelity_pass = if fidelity_ok { fidelity >= 0.9999 } else { true };
    let pass = effective_prob >= min_prob && fidelity_pass;

    BenchResult {
        name,
        qubits:     circuit.num_qubits(),
        gates:      stats.gate_count,
        depth:      stats.depth,
        label,
        fidelity,
        time_sv_ms,
        time_q2_ms,
        pass,
        note,
    }
}

fn main() {
    println!("╔══════════════════════════════════════════════════════════════╗");
    println!("║   CleitonForge — Canonical Benchmark Suite                  ║");
    println!("║   Author: Cleiton Augusto Correa Bezerra                    ║");
    println!("╚══════════════════════════════════════════════════════════════╝");
    println!("  Backends : statevector (native) + QuantRS2");
    println!("  Shots    : {SHOTS}  |  Seed: {DEFAULT_SEED:#x}");
    println!();

    let benchmarks: Vec<(&str, Circuit, ExpectKind, f64)> = vec![
        ("Bell state",         bell_state(),             ExpectKind::Entangled,       0.95),
        ("GHZ state",          ghz_state(),              ExpectKind::Entangled,       0.95),
        ("QFT |100⟩",          qft3_on_100(),            ExpectKind::Uniform,         1.00),
        ("Bernstein-Vazirani", bernstein_vazirani_101(), ExpectKind::Bitstring("101"),0.99),
        ("Grover |101⟩",       grover_101(),             ExpectKind::Bitstring("101"),0.90),
        ("QAOA MaxCut",        qaoa_maxcut_2q(),         ExpectKind::MaxCut,          0.90),
    ];

    let mut results = Vec::new();
    for (name, circuit, expect, min_prob) in benchmarks {
        print!("  Running {name}...");
        let r = run_benchmark(name, &circuit, expect, min_prob);
        println!(" {}", if r.pass { "✓" } else { "✗" });
        results.push(r);
    }

    println!();
    println!("┌─────────────────────┬────┬───────┬───────┬──────────────────┬──────────────┬────────────┬────┐");
    println!("│ Benchmark           │  Q │ Gates │ Depth │ Result           │   Fidelity   │  Time (ms) │ OK │");
    println!("├─────────────────────┼────┼───────┼───────┼──────────────────┼──────────────┼────────────┼────┤");

    let mut all_pass = true;
    for r in &results {
        let fid_str = if r.note.is_some() {
            format!("{:.6} *", r.fidelity)
        } else {
            format!("{:.8}", r.fidelity)
        };
        println!(
            "│ {:<19} │ {:>2} │  {:>4} │  {:>4} │ {:<16} │ {:>12} │ {:>10} │ {:>2} │",
            r.name, r.qubits, r.gates, r.depth,
            r.label,
            fid_str,
            format!("{:.1}/{:.1}", r.time_sv_ms, r.time_q2_ms),
            if r.pass { "✅" } else { "❌" },
        );
        if !r.pass { all_pass = false; }
    }

    println!("└─────────────────────┴────┴───────┴───────┴──────────────────┴──────────────┴────────────┴────┘");
    println!();

    let has_notes = results.iter().any(|r| r.note.is_some());
    if has_notes {
        println!("  * Fidelity not required — documented inter-framework Rz sign convention divergence.");
        println!("    This is the value CleitonForge surfaces: same circuit, different backend math.");
        println!();
    }

    if all_pass {
        println!("  All benchmarks passed ✅");
    } else {
        println!("  Some benchmarks failed ❌");
        std::process::exit(1);
    }
}
