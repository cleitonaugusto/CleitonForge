//! `cforge bench` — the canonical benchmark suite as a CI gate.
//!
//! Runs six standard algorithms across the selected backends, checks each one
//! against a known expected outcome, and compares the backends against each
//! other. Exits 1 if anything fails, so it can be dropped straight into a
//! pipeline step.
//!
//! Cross-backend fidelity is only meaningful with two or more backends. With a
//! single backend the suite still checks every expected outcome; the fidelity
//! column simply reports nothing to compare.

use std::time::Instant;

use comfy_table::{presets::UTF8_FULL, Table};
use serde_json::{json, Map, Value};

use cforge_backends::SimulationBackend;
use cforge_core::{Circuit, GateKind, Operation};
use cforge_metrics::{compute_stats, statevector_fidelity};

// ── Circuit builders ──────────────────────────────────────────────────────────

fn bell_state() -> Circuit {
    let mut c = Circuit::new(2);
    c.push(Operation::new(GateKind::H, vec![0], vec![]));
    c.push(Operation::new(GateKind::Cx, vec![0, 1], vec![]));
    c
}

fn ghz_state() -> Circuit {
    let mut c = Circuit::new(3);
    c.push(Operation::new(GateKind::H, vec![0], vec![]));
    c.push(Operation::new(GateKind::Cx, vec![0, 1], vec![]));
    c.push(Operation::new(GateKind::Cx, vec![1, 2], vec![]));
    c
}

fn qft3_on_100() -> Circuit {
    use std::f64::consts::PI;
    let mut c = Circuit::new(3);
    // Initialise to |100⟩
    c.push(Operation::new(GateKind::X, vec![0], vec![]));
    // QFT
    c.push(Operation::new(GateKind::H, vec![0], vec![]));
    c.push(Operation::new(GateKind::Cp, vec![1, 0], vec![PI / 2.0]));
    c.push(Operation::new(GateKind::Cp, vec![2, 0], vec![PI / 4.0]));
    c.push(Operation::new(GateKind::H, vec![1], vec![]));
    c.push(Operation::new(GateKind::Cp, vec![2, 1], vec![PI / 2.0]));
    c.push(Operation::new(GateKind::H, vec![2], vec![]));
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
        c.push(Operation::new(GateKind::X, vec![1], vec![]));
        c.push(Operation::new(GateKind::H, vec![2], vec![]));
        c.push(Operation::new(GateKind::Ccx, vec![0, 1, 2], vec![]));
        c.push(Operation::new(GateKind::H, vec![2], vec![]));
        c.push(Operation::new(GateKind::X, vec![1], vec![]));
        // Diffuser
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

// QAOA p=1 MaxCut on a single edge (q0-q1).
// Optimal angles (statevector): γ = -3π/4, β = -π/8 → 100% |01⟩+|10⟩ (cut states).
// Note: quantrs2 applies Rz with the opposite imaginary-part sign convention,
// so the two backends intentionally diverge here. This is one of the inter-
// framework discrepancies CleitonForge is designed to surface.
fn qaoa_maxcut_2q() -> Circuit {
    use std::f64::consts::PI;
    let gamma = -3.0 * PI / 4.0;
    let beta = -PI / 8.0;
    let mut c = Circuit::new(2);
    c.push(Operation::new(GateKind::H, vec![0], vec![]));
    c.push(Operation::new(GateKind::H, vec![1], vec![]));
    // Cost: RZZ(2γ) = CX, RZ(2γ), CX
    c.push(Operation::new(GateKind::Cx, vec![0, 1], vec![]));
    c.push(Operation::new(GateKind::Rz, vec![1], vec![2.0 * gamma]));
    c.push(Operation::new(GateKind::Cx, vec![0, 1], vec![]));
    // Mixer: RX(2β)
    c.push(Operation::new(GateKind::Rx, vec![0], vec![2.0 * beta]));
    c.push(Operation::new(GateKind::Rx, vec![1], vec![2.0 * beta]));
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

/// Knobs that apply to every case in a run.
struct BenchOptions {
    shots: usize,
    seed: u64,
    threshold: f64,
}

#[derive(Debug)]
struct BenchResult {
    name: &'static str,
    qubits: usize,
    gates: usize,
    depth: usize,
    label: String,
    outcome_prob: f64,
    min_prob: f64,
    /// `None` when there is nothing to compare against: either a single backend
    /// was selected, or the case documents a deliberate divergence.
    fidelity: Option<f64>,
    times_ms: Vec<(String, f64)>,
    pass: bool,
    note: Option<&'static str>,
}

/// Runs one case across every selected backend.
///
/// The first backend is the reference: expected outcomes are checked against
/// its counts, and every other backend is compared to it by state fidelity.
fn run_benchmark(
    name: &'static str,
    circuit: &Circuit,
    expect: ExpectKind,
    min_prob: f64,
    backends: &[Box<dyn SimulationBackend>],
    opts: &BenchOptions,
) -> Result<BenchResult, String> {
    let stats = compute_stats(circuit);

    let mut times_ms = Vec::with_capacity(backends.len());
    let mut runs = Vec::with_capacity(backends.len());
    for b in backends {
        let t0 = Instant::now();
        let out = b
            .run(circuit, opts.shots, opts.seed)
            .map_err(|e| format!("{name}: backend {} failed: {e}", b.name()))?;
        times_ms.push((b.name().to_string(), t0.elapsed().as_secs_f64() * 1000.0));
        runs.push(out);
    }

    let reference = &runs[0];
    let n = circuit.num_qubits();
    let total = reference.counts.values().sum::<usize>() as f64;
    if total == 0.0 {
        // Every expected outcome here is a ratio over sampled counts. With no
        // shots the denominator is zero and each ratio comes back NaN, which
        // compares false against its threshold and fails the case for a reason
        // that has nothing to do with the code under test.
        return Err(format!(
            "{name}: backend {} returned no shots; bench needs --shots >= 1",
            backends[0].name()
        ));
    }
    let mut note = None;

    let (outcome_prob, label, compare_backends) = match &expect {
        ExpectKind::Entangled => {
            let zeros = "0".repeat(n);
            let ones = "1".repeat(n);
            let p0 = reference.counts.get(zeros.as_str()).copied().unwrap_or(0) as f64 / total;
            let p1 = reference.counts.get(ones.as_str()).copied().unwrap_or(0) as f64 / total;
            (
                p0 + p1,
                format!("|{zeros}⟩+|{ones}⟩ {:.0}%", (p0 + p1) * 100.0),
                true,
            )
        }
        ExpectKind::Bitstring(s) => {
            let p = reference.counts.get(*s).copied().unwrap_or(0) as f64 / total;
            (p, format!("|{s}⟩ {:.0}%", p * 100.0), true)
        }
        ExpectKind::Uniform => {
            let expected_p = 1.0 / (1 << n) as f64;
            let uniform = reference
                .statevector
                .iter()
                .all(|a| (a.norm_sqr() - expected_p).abs() < 0.002);
            let lbl = if uniform {
                format!("uniform 1/{} ✓", 1 << n)
            } else {
                "not uniform".to_string()
            };
            (if uniform { 1.0 } else { 0.0 }, lbl, true)
        }
        ExpectKind::MaxCut => {
            let p01 = reference.counts.get("01").copied().unwrap_or(0) as f64 / total;
            let p10 = reference.counts.get("10").copied().unwrap_or(0) as f64 / total;
            let cut = p01 + p10;
            note = Some("Rz convention divergence (by design)");
            // Which way this lands depends on the reference backend's Rz sign
            // convention: the native statevector concentrates on the cut states,
            // quantrs2 on their complement. Both are correct for their own
            // convention, so the check is that the outcome is one of the two —
            // a result in the middle is a real regression in either.
            let decisive = cut.max(1.0 - cut);
            (
                decisive,
                format!("cut {:.0}% ({})", cut * 100.0, if cut >= 0.5 { "cut states" } else { "complement" }),
                false,
            )
        }
    };

    // Worst fidelity of any backend against the reference. With one backend
    // there is no pair and so no number to report. It is always computed when
    // there is one to compute, including for the deliberate-divergence case:
    // that value is the disagreement itself, which is worth seeing even though
    // it must not fail the build.
    let fidelity = if runs.len() > 1 {
        Some(
            runs[1..]
                .iter()
                .map(|r| {
                    statevector_fidelity(&reference.statevector, &r.statevector).unwrap_or(0.0)
                })
                .fold(f64::INFINITY, f64::min),
        )
    } else {
        None
    };

    let fidelity_pass = !compare_backends || fidelity.is_none_or(|f| f >= opts.threshold);
    let pass = outcome_prob >= min_prob && fidelity_pass;

    Ok(BenchResult {
        name,
        qubits: circuit.num_qubits(),
        gates: stats.gate_count,
        depth: stats.depth,
        label,
        outcome_prob,
        min_prob,
        fidelity,
        times_ms,
        pass,
        note,
    })
}

/// The suite itself: name, circuit, expected outcome, and the minimum
/// probability that outcome must reach before the case counts as passing.
fn suite() -> Vec<(&'static str, Circuit, ExpectKind, f64)> {
    vec![
        ("Bell state", bell_state(), ExpectKind::Entangled, 0.95),
        ("GHZ state", ghz_state(), ExpectKind::Entangled, 0.95),
        ("QFT |100⟩", qft3_on_100(), ExpectKind::Uniform, 1.00),
        (
            "Bernstein-Vazirani",
            bernstein_vazirani_101(),
            ExpectKind::Bitstring("101"),
            0.99,
        ),
        (
            "Grover |101⟩",
            grover_101(),
            ExpectKind::Bitstring("101"),
            0.90,
        ),
        ("QAOA MaxCut", qaoa_maxcut_2q(), ExpectKind::MaxCut, 0.90),
    ]
}

fn to_json(results: &[BenchResult], names: &[String], shots: usize, seed: u64, threshold: f64) -> Value {
    let cases: Vec<Value> = results
        .iter()
        .map(|r| {
            let mut times = Map::new();
            for (backend, ms) in &r.times_ms {
                times.insert(backend.clone(), json!(ms));
            }
            json!({
                "name": r.name,
                "qubits": r.qubits,
                "gates": r.gates,
                "depth": r.depth,
                "outcome": r.label,
                "outcome_probability": r.outcome_prob,
                "min_probability": r.min_prob,
                "fidelity": r.fidelity,
                "time_ms": times,
                "passed": r.pass,
                "note": r.note,
            })
        })
        .collect();

    let failed = results.iter().filter(|r| !r.pass).count();
    // Only gated cases belong in this number. Including the documented Rz
    // divergence would pin worst_fidelity at ~0 on every green run, and a field
    // that reads the same whether or not anything is wrong is one a dashboard
    // learns to ignore.
    let worst = results
        .iter()
        .filter(|r| r.note.is_none())
        .filter_map(|r| r.fidelity.map(|f| (r.name, f)))
        .min_by(|a, b| a.1.total_cmp(&b.1));
    let documented: Vec<Value> = results
        .iter()
        .filter(|r| r.note.is_some())
        .map(|r| json!({"name": r.name, "fidelity": r.fidelity, "note": r.note}))
        .collect();

    json!({
        "backends": names,
        "seed": seed,
        "shots": shots,
        "threshold": threshold,
        "results": cases,
        "summary": {
            "total": results.len(),
            "passed": results.len() - failed,
            "failed": failed,
            "worst_fidelity": worst.map(|(_, f)| f),
            "worst_circuit": worst.map(|(n, _)| n),
            "documented_divergences": documented,
        }
    })
}

fn print_table(results: &[BenchResult], names: &[String], shots: usize, seed: u64, threshold: f64) {
    println!("╔══════════════════════════════════════════════════════════════╗");
    println!("║        CleitonForge — Canonical Benchmark Suite              ║");
    println!("╠══════════════════════════════════════════════════════════════╣");
    println!("║  Backends : {:<49}║", names.join(", "));
    println!("║  Shots    : {:<49}║", format!("{shots}  |  seed {seed:#x}"));
    println!("║  Threshold: {:<49}║", format!("fidelity ≥ {threshold}"));
    println!("╚══════════════════════════════════════════════════════════════╝");
    println!();

    let mut table = Table::new();
    table.load_preset(UTF8_FULL);
    table.set_header(vec![
        "Benchmark", "Q", "Gates", "Depth", "Result", "Fidelity", "Time (ms)", "OK",
    ]);

    for r in results {
        let fidelity = match (r.fidelity, r.note) {
            (Some(f), Some(_)) => format!("{f:.6} *"),
            (Some(f), None) => format!("{f:.8}"),
            (None, _) => "—".to_string(),
        };
        let times = r
            .times_ms
            .iter()
            .map(|(_, ms)| format!("{ms:.1}"))
            .collect::<Vec<_>>()
            .join("/");
        table.add_row(vec![
            r.name.to_string(),
            r.qubits.to_string(),
            r.gates.to_string(),
            r.depth.to_string(),
            r.label.clone(),
            fidelity,
            times,
            if r.pass { "✅" } else { "❌" }.to_string(),
        ]);
    }
    println!("{table}");

    // Only worth explaining the asterisk if a row actually carries one.
    if results.iter().any(|r| r.note.is_some() && r.fidelity.is_some()) {
        println!();
        println!("  * Fidelity not gated — documented inter-framework Rz sign convention");
        println!("    divergence. Same circuit, different backend math: that disagreement");
        println!("    is the thing CleitonForge exists to surface, not a regression.");
    }
    if results.iter().any(|r| r.fidelity.is_none() && r.note.is_none()) {
        println!();
        println!("  — No cross-backend fidelity: only one backend was selected. Expected");
        println!("    outcomes were still checked. Pass two or more to gate on agreement.");
    }
}

/// Entry point for `cforge bench`. Returns the process exit code.
pub fn cmd_bench(
    backends: Vec<Box<dyn SimulationBackend>>,
    shots: usize,
    seed: u64,
    threshold: f64,
    json_output: bool,
) -> i32 {
    let names: Vec<String> = backends.iter().map(|b| b.name().to_string()).collect();
    let opts = BenchOptions {
        shots,
        seed,
        threshold,
    };

    let mut results = Vec::new();
    for (name, circuit, expect, min_prob) in suite() {
        if !json_output {
            print!("  Running {name}...");
        }
        match run_benchmark(name, &circuit, expect, min_prob, &backends, &opts) {
            Ok(r) => {
                if !json_output {
                    println!(" {}", if r.pass { "✓" } else { "✗" });
                }
                results.push(r);
            }
            Err(e) => {
                eprintln!("error: {e}");
                return 2;
            }
        }
    }

    if json_output {
        println!(
            "{}",
            serde_json::to_string_pretty(&to_json(&results, &names, shots, seed, threshold))
                .expect("serialising the report cannot fail")
        );
    } else {
        println!();
        print_table(&results, &names, shots, seed, threshold);
        println!();
    }

    let failed = results.iter().filter(|r| !r.pass).count();
    if failed == 0 {
        if !json_output {
            println!("  All {} benchmarks passed ✅", results.len());
        }
        0
    } else {
        if !json_output {
            println!("  {failed} of {} benchmarks failed ❌", results.len());
        }
        1
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cforge_backends::{NativeStateVectorBackend, QuantRS2Backend, DEFAULT_SEED};

    fn native() -> Vec<Box<dyn SimulationBackend>> {
        vec![Box::new(NativeStateVectorBackend)]
    }

    fn both() -> Vec<Box<dyn SimulationBackend>> {
        vec![
            Box::new(NativeStateVectorBackend),
            Box::new(QuantRS2Backend),
        ]
    }

    fn run_one(
        name: &'static str,
        backends: &[Box<dyn SimulationBackend>],
        threshold: f64,
    ) -> BenchResult {
        let (n, circuit, expect, min_prob) = suite()
            .into_iter()
            .find(|(n, ..)| *n == name)
            .expect("benchmark is in the suite");
        let opts = BenchOptions {
            shots: 4096,
            seed: DEFAULT_SEED,
            threshold,
        };
        run_benchmark(n, &circuit, expect, min_prob, backends, &opts).expect("backends run")
    }

    #[test]
    fn suite_has_six_benchmarks() {
        assert_eq!(suite().len(), 6);
    }

    #[test]
    fn every_benchmark_passes_on_both_backends() {
        for (name, ..) in suite() {
            let r = run_one(name, &both(), 0.9999);
            assert!(r.pass, "{name} should pass with the default threshold");
        }
    }

    /// A single backend has nothing to compare against, so no fidelity is
    /// reported — but the expected outcome is still checked.
    #[test]
    fn single_backend_reports_no_fidelity_but_still_checks_outcome() {
        let r = run_one("Bell state", &native(), 0.9999);
        assert!(r.fidelity.is_none());
        assert!(r.pass);
        assert!(r.outcome_prob >= r.min_prob);
    }

    /// The gate has to be able to fail. A threshold no fidelity can reach must
    /// turn a compared benchmark red, or the whole CI step is decoration.
    #[test]
    fn unreachable_threshold_fails_a_compared_benchmark() {
        let r = run_one("Bell state", &both(), 1.1);
        assert!(!r.pass, "fidelity 1.0 must not clear a threshold of 1.1");
    }

    /// ...and it must fail for the right reason: with one backend there is no
    /// fidelity, so the same impossible threshold must not fail anything.
    #[test]
    fn unreachable_threshold_is_inert_without_a_second_backend() {
        let r = run_one("Bell state", &native(), 1.1);
        assert!(r.pass);
    }

    fn reversed() -> Vec<Box<dyn SimulationBackend>> {
        vec![
            Box::new(QuantRS2Backend),
            Box::new(NativeStateVectorBackend),
        ]
    }

    /// Which backend the user happens to list first must not decide whether the
    /// build is red. The QAOA case lands on the cut states under one Rz
    /// convention and on their complement under the other; both are correct.
    #[test]
    fn backend_order_does_not_change_the_verdict() {
        for (name, ..) in suite() {
            assert!(
                run_one(name, &reversed(), 0.9999).pass,
                "{name} must pass with quantrs2 as the reference backend"
            );
        }
    }

    /// Same argument with no second backend at all: a CI box that only has
    /// quantrs2 installed must not see a permanently failing suite.
    #[test]
    fn quantrs2_alone_passes_every_benchmark() {
        let only = vec![Box::new(QuantRS2Backend) as Box<dyn SimulationBackend>];
        for (name, ..) in suite() {
            assert!(run_one(name, &only, 0.9999).pass, "{name} on quantrs2 alone");
        }
    }

    /// Zero shots leaves the counts empty, so every expected outcome would be
    /// 0/0. That must be a loud error, not a NaN that quietly compares false
    /// and reddens the build for the wrong reason.
    #[test]
    fn zero_shots_is_an_error_not_a_silent_nan() {
        let backends = native();
        let (n, circuit, expect, min_prob) = suite().into_iter().next().unwrap();
        let opts = BenchOptions {
            shots: 0,
            seed: DEFAULT_SEED,
            threshold: 0.9999,
        };
        let err = run_benchmark(n, &circuit, expect, min_prob, &backends, &opts)
            .expect_err("zero shots must not produce a result");
        assert!(err.contains("--shots"), "the error should say how to fix it");
    }

    /// `worst_fidelity` gates dashboards, so it must describe the cases that are
    /// actually gated. Including the documented divergence would pin it near
    /// zero on every green run.
    #[test]
    fn summary_worst_fidelity_ignores_the_documented_divergence() {
        let backends = both();
        let opts = BenchOptions {
            shots: 4096,
            seed: DEFAULT_SEED,
            threshold: 0.9999,
        };
        let results: Vec<BenchResult> = suite()
            .into_iter()
            .map(|(n, c, e, m)| run_benchmark(n, &c, e, m, &backends, &opts).unwrap())
            .collect();
        let names = vec!["a".to_string(), "b".to_string()];
        let v = to_json(&results, &names, 4096, DEFAULT_SEED, 0.9999);
        let worst = v["summary"]["worst_fidelity"].as_f64().unwrap();
        assert!(worst > 0.99, "worst gated fidelity was {worst}");
        assert_eq!(
            v["summary"]["documented_divergences"]
                .as_array()
                .unwrap()
                .len(),
            1,
            "the divergence must still be reported, just not as the worst case"
        );
    }

    /// The Rz convention divergence is documented, not a regression: its
    /// fidelity is reported but must never gate the build.
    #[test]
    fn maxcut_divergence_is_reported_but_not_gated() {
        let r = run_one("QAOA MaxCut", &both(), 1.1);
        assert!(r.note.is_some(), "the divergence must be annotated");
        assert!(r.fidelity.is_some(), "the divergence value is worth seeing");
        assert!(r.pass, "a documented divergence must not fail the build");
    }
}
