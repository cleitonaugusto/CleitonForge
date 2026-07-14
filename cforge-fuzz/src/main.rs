use std::f64::consts::PI;

use clap::Parser;
use num_complex::Complex64;
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;
use rayon::prelude::*;

use cforge_backends::{NativeStateVectorBackend, QuantRS2Backend, SimulationBackend};
use cforge_core::{Circuit, GateKind, Operation};

// ── CLI ──────────────────────────────────────────────────────────────────────

#[derive(Parser)]
#[command(
    name = "cforge-fuzz",
    about = "Quantum differential fuzzer — finds gate convention divergences automatically"
)]
struct Cli {
    /// Number of random circuits to generate and test
    #[arg(long, default_value_t = 10_000)]
    iterations: usize,

    /// Maximum qubits per generated circuit
    #[arg(long, default_value_t = 3)]
    max_qubits: usize,

    /// Maximum gates per generated circuit
    #[arg(long, default_value_t = 8)]
    max_depth: usize,

    /// PRNG seed (0 = random)
    #[arg(long, default_value_t = 0)]
    seed: u64,

    /// Amplitude tolerance for divergence detection
    #[arg(long, default_value_t = 1e-6)]
    tol: f64,

    /// Stop after finding this many divergences (0 = keep going)
    #[arg(long, default_value_t = 1)]
    stop_after: usize,
}

// ── Random circuit generation ─────────────────────────────────────────────────

/// Gates that expose Rz sign and phase convention bugs.
/// Deliberately excludes Clifford-only gates that are convention-safe.
const CONVENTION_SENSITIVE: &[GateKind] = &[
    GateKind::H,
    GateKind::Rz,
    GateKind::Rx,
    GateKind::Ry,
    GateKind::Phase,
    GateKind::T,
    GateKind::Tdg,
    GateKind::S,
    GateKind::Sdg,
    GateKind::Cx,
    GateKind::Cz,
    GateKind::Cp,
    GateKind::Crz,
];

fn random_circuit(rng: &mut StdRng, num_qubits: usize, depth: usize) -> Circuit {
    let mut c = Circuit::new(num_qubits);
    for _ in 0..depth {
        let gate = CONVENTION_SENSITIVE[rng.gen_range(0..CONVENTION_SENSITIVE.len())];
        let nq = gate.num_qubits().min(num_qubits);
        if nq == 0 || num_qubits < nq { continue; }

        // Pick distinct qubits
        let qubits: Vec<usize> = {
            let mut available: Vec<usize> = (0..num_qubits).collect();
            let mut picked = Vec::with_capacity(nq);
            for _ in 0..nq {
                let idx = rng.gen_range(0..available.len());
                picked.push(available.remove(idx));
            }
            picked
        };

        let params: Vec<f64> = (0..gate.num_params())
            .map(|_| rng.gen_range(0.0..2.0 * PI))
            .collect();

        c.push(Operation::new(gate, qubits, params));
    }
    c
}

// ── Differential oracle ───────────────────────────────────────────────────────

struct Divergence {
    circuit:   Circuit,
    sv_a:      Vec<Complex64>,
    sv_b:      Vec<Complex64>,
    max_delta: f64,
}

fn amplitude_distance(a: &[Complex64], b: &[Complex64]) -> f64 {
    if a.len() != b.len() { return f64::INFINITY; }
    a.iter().zip(b).map(|(x, y)| (x - y).norm()).fold(0.0_f64, f64::max)
}

fn run_backend(backend: &dyn SimulationBackend, circuit: &Circuit) -> Option<Vec<Complex64>> {
    backend.run(circuit, 0, 0xdeadbeef).ok().map(|r| r.statevector)
}

fn check_divergence(
    circuit: &Circuit,
    tol: f64,
) -> Option<Divergence> {
    let sv_a = run_backend(&NativeStateVectorBackend, circuit)?;
    let sv_b = run_backend(&QuantRS2Backend, circuit)?;
    let delta = amplitude_distance(&sv_a, &sv_b);
    if delta > tol {
        Some(Divergence { circuit: circuit.clone(), sv_a, sv_b, max_delta: delta })
    } else {
        None
    }
}

// ── Delta-debugging minimizer ─────────────────────────────────────────────────

/// Reduce the circuit to the smallest subsequence that still diverges.
fn minimize(original: &Circuit, tol: f64) -> Circuit {
    let ops: Vec<Operation> = original.operations.iter().cloned().collect();
    let n = ops.len();
    let nq = original.num_qubits();

    // Try removing each gate one at a time (greedy single-pass)
    let mut current = ops.clone();
    let mut changed = true;
    while changed {
        changed = false;
        let mut i = 0;
        while i < current.len() {
            let mut candidate = current.clone();
            candidate.remove(i);
            let mut c = Circuit::new(nq);
            for op in &candidate { c.push(op.clone()); }
            if check_divergence(&c, tol).is_some() {
                current = candidate;
                changed = true;
                // don't increment i — next element shifted down
            } else {
                i += 1;
            }
        }
    }

    let mut minimal = Circuit::new(nq);
    for op in &current { minimal.push(op.clone()); }
    minimal
}

// ── Reporting ─────────────────────────────────────────────────────────────────

fn print_divergence(d: &Divergence, minimal: &Circuit) {
    println!();
    println!("╔══════════════════════════════════════════════════════════════╗");
    println!("║                 DIVERGENCE FOUND                            ║");
    println!("╚══════════════════════════════════════════════════════════════╝");
    println!();
    println!("Original circuit : {} gates, {} qubits", d.circuit.operations.len(), d.circuit.num_qubits());
    println!("Minimal witness  : {} gates", minimal.operations.len());
    println!("Max |Δamplitude| : {:.2e}", d.max_delta);
    println!();
    println!("Minimal counterexample:");
    for op in &minimal.operations {
        let params_str = if op.params.is_empty() {
            String::new()
        } else {
            format!("({})", op.params.iter().map(|p| format!("{:.4}", p)).collect::<Vec<_>>().join(", "))
        };
        println!("  {} q[{}]{}", op.kind.qasm_name(), op.qubits.iter().map(|q| q.to_string()).collect::<Vec<_>>().join(",q["), params_str);
    }
    println!();

    println!("Statevector comparison (first 8 amplitudes):");
    println!("  {:>5}  {:>30}  {:>30}", "index", "statevector (native)", "statevector (quantrs2)");
    let len = d.sv_a.len().min(8);
    for i in 0..len {
        let a = d.sv_a[i];
        let b = d.sv_b[i];
        let delta = (a - b).norm();
        let flag = if delta > 1e-9 { " ←" } else { "" };
        println!("  [{:>3}]  {:>+.6}+{:>+.6}i  {:>+.6}+{:>+.6}i{}", i, a.re, a.im, b.re, b.im, flag);
    }
    println!();
    println!("Likely cause: Rz sign inversion (native: diag(e^{{-iλ/2}}, e^{{+iλ/2}}); quantrs2: reversed)");
    println!();
}

// ── Main ──────────────────────────────────────────────────────────────────────

fn main() {
    let cli = Cli::parse();

    let seed = if cli.seed == 0 {
        rand::thread_rng().gen()
    } else {
        cli.seed
    };

    println!("╔══════════════════════════════════════════════════════════════╗");
    println!("║    cforge-fuzz — Quantum Differential Fuzzer                ║");
    println!("╠══════════════════════════════════════════════════════════════╣");
    println!("║  Backends : statevector (oracle) vs quantrs2                ║");
    println!("║  Seed     : {:<48}║", format!("0x{:x}", seed));
    println!("║  Circuits : {:<48}║", cli.iterations);
    println!("║  Tolerance: {:<48}║", format!("{:.0e}", cli.tol));
    println!("╚══════════════════════════════════════════════════════════════╝");
    println!();

    let mut found = 0usize;

    for i in 0..cli.iterations {
        if i % 1000 == 0 && i > 0 {
            eprintln!("[{}/{}] {} divergences so far...", i, cli.iterations, found);
        }

        let mut rng = StdRng::seed_from_u64(seed.wrapping_add(i as u64));
        let nq = rng.gen_range(1..=cli.max_qubits);
        let depth = rng.gen_range(1..=cli.max_depth);
        let circuit = random_circuit(&mut rng, nq, depth);

        if let Some(div) = check_divergence(&circuit, cli.tol) {
            found += 1;
            let minimal = minimize(&div.circuit, cli.tol);
            print_divergence(&div, &minimal);

            if cli.stop_after > 0 && found >= cli.stop_after {
                println!("Stopping after {} divergence(s). Use --stop-after 0 to continue.", found);
                break;
            }
        }
    }

    println!();
    if found == 0 {
        println!("Result: No divergences found in {} circuits. ✅", cli.iterations);
    } else {
        println!("Result: {} divergence(s) found in {} circuits. ❌", found, cli.iterations);
        std::process::exit(1);
    }
}
