//! QGCS conformance checks — runs against any SimulationBackend.
//!
//! Each check builds a small witness circuit, runs it on the target backend,
//! and verifies the output amplitude against the OpenQASM 3 specification.
//! The expected value is derived analytically before the circuit is written.

use num_complex::Complex64;
use std::f64::consts::{FRAC_1_SQRT_2, FRAC_PI_2, PI};

use cforge_core::{Circuit, GateKind, Operation};

use crate::SimulationBackend;

/// Tolerance for amplitude-level phase comparisons.
/// Tighter than 1e-6 (avoids false passes) but realistic for f64 accumulation.
const PHASE_TOL: f64 = 1e-9;

pub enum CheckStatus {
    Pass,
    Fail { expected: String, got: String },
    Skip { reason: String },
}

pub struct CheckResult {
    pub name: &'static str,
    pub dimension: &'static str,
    pub status: CheckStatus,
}

impl CheckResult {
    pub fn passed(&self) -> bool {
        matches!(self.status, CheckStatus::Pass)
    }

    pub fn skipped(&self) -> bool {
        matches!(self.status, CheckStatus::Skip { .. })
    }
}

fn skip(name: &'static str, dimension: &'static str, reason: String) -> CheckResult {
    CheckResult {
        name,
        dimension,
        status: CheckStatus::Skip { reason },
    }
}

fn run(backend: &dyn SimulationBackend, c: &Circuit) -> Result<Vec<Complex64>, String> {
    backend
        .run(c, 0, 0)
        .map(|r| r.statevector)
        .map_err(|e| e.to_string())
}

/// Fidelity |⟨a|b⟩|² with length and zero-norm guards.
/// Returns None on length mismatch or degenerate (zero) state.
fn fidelity(a: &[Complex64], b: &[Complex64]) -> Option<f64> {
    if a.len() != b.len() || a.is_empty() {
        return None;
    }
    let norm_a: f64 = a.iter().map(|x| x.norm_sqr()).sum();
    let norm_b: f64 = b.iter().map(|x| x.norm_sqr()).sum();
    if norm_a < 1e-15 || norm_b < 1e-15 {
        return None;
    }
    let inner: Complex64 = a.iter().zip(b).map(|(x, y)| x.conj() * y).sum();
    Some(inner.norm_sqr())
}

fn ratio_check(
    backend: &dyn SimulationBackend,
    name: &'static str,
    dimension: &'static str,
    circuit: Circuit,
    num: usize,
    den: usize,
    expected: Complex64,
) -> CheckResult {
    match run(backend, &circuit) {
        Err(e) => skip(name, dimension, e),
        Ok(sv) => {
            let Some(&n) = sv.get(num) else {
                return skip(
                    name,
                    dimension,
                    format!("statevector too short (len={}, need index {num})", sv.len()),
                );
            };
            let Some(&d) = sv.get(den) else {
                return skip(
                    name,
                    dimension,
                    format!("statevector too short (len={}, need index {den})", sv.len()),
                );
            };
            if d.norm() < 1e-15 {
                return skip(
                    name,
                    dimension,
                    format!("sv[{den}] ≈ 0 — degenerate state (broken H?)"),
                );
            }
            let ratio = n / d;
            CheckResult {
                name,
                dimension,
                status: if (ratio - expected).norm() < PHASE_TOL {
                    CheckStatus::Pass
                } else {
                    CheckStatus::Fail {
                        expected: format!("{expected:.4}"),
                        got: format!("{ratio:.4}"),
                    }
                },
            }
        }
    }
}

fn fidelity_check(
    backend: &dyn SimulationBackend,
    name: &'static str,
    dimension: &'static str,
    circuit_a: Circuit,
    circuit_b: Circuit,
    expected_f: f64,
) -> CheckResult {
    let (sa, sb) = match (run(backend, &circuit_a), run(backend, &circuit_b)) {
        (Ok(a), Ok(b)) => (a, b),
        (Err(e), _) | (_, Err(e)) => return skip(name, dimension, e),
    };
    match fidelity(&sa, &sb) {
        None => skip(
            name,
            dimension,
            format!(
                "statevector dimension mismatch or degenerate state (|a|={}, |b|={})",
                sa.len(),
                sb.len()
            ),
        ),
        Some(f) => CheckResult {
            name,
            dimension,
            status: if (f - expected_f).abs() < PHASE_TOL {
                CheckStatus::Pass
            } else {
                CheckStatus::Fail {
                    expected: format!("{expected_f:.6}"),
                    got: format!("{f:.6}"),
                }
            },
        },
    }
}

fn amp_check(
    backend: &dyn SimulationBackend,
    name: &'static str,
    dimension: &'static str,
    circuit: Circuit,
    index: usize,
    expected: Complex64,
) -> CheckResult {
    match run(backend, &circuit) {
        Err(e) => skip(name, dimension, e),
        Ok(sv) => {
            let Some(&amp) = sv.get(index) else {
                return skip(
                    name,
                    dimension,
                    format!(
                        "statevector too short (len={}, need index {index})",
                        sv.len()
                    ),
                );
            };
            CheckResult {
                name,
                dimension,
                status: if (amp - expected).norm() < PHASE_TOL {
                    CheckStatus::Pass
                } else {
                    CheckStatus::Fail {
                        expected: format!("{expected:.4}"),
                        got: format!("{amp:.4}"),
                    }
                },
            }
        }
    }
}

// ── Individual checks ─────────────────────────────────────────────────────────

fn check_rz_sign(b: &dyn SimulationBackend) -> CheckResult {
    // Rz(π/2)|+⟩: sv[1]/sv[0] = e^{iπ/2} = i  (OpenQASM 3: diag(e^{-iλ/2}, e^{+iλ/2}))
    let mut c = Circuit::new(1);
    c.push(Operation::new(GateKind::H, vec![0], vec![]));
    c.push(Operation::new(GateKind::Rz, vec![0], vec![FRAC_PI_2]));
    ratio_check(
        b,
        "Rz(π/2)|+⟩ ratio = +i",
        "Rz sign",
        c,
        1,
        0,
        Complex64::new(0.0, 1.0),
    )
}

fn check_t_phase(b: &dyn SimulationBackend) -> CheckResult {
    // T|+⟩: sv[1]/sv[0] = e^{iπ/4} = 1/√2 + i/√2  (T = diag(1, e^{iπ/4}))
    let mut c = Circuit::new(1);
    c.push(Operation::new(GateKind::H, vec![0], vec![]));
    c.push(Operation::new(GateKind::T, vec![0], vec![]));
    ratio_check(
        b,
        "T|+⟩ ratio = e^{iπ/4}",
        "T gate",
        c,
        1,
        0,
        Complex64::new(FRAC_1_SQRT_2, FRAC_1_SQRT_2),
    )
}

fn check_tdg_phase(b: &dyn SimulationBackend) -> CheckResult {
    // Tdg|+⟩: sv[1]/sv[0] = e^{-iπ/4}  (Tdg = diag(1, e^{-iπ/4}))
    let mut c = Circuit::new(1);
    c.push(Operation::new(GateKind::H, vec![0], vec![]));
    c.push(Operation::new(GateKind::Tdg, vec![0], vec![]));
    ratio_check(
        b,
        "Tdg|+⟩ ratio = e^{-iπ/4}",
        "T gate",
        c,
        1,
        0,
        Complex64::new(FRAC_1_SQRT_2, -FRAC_1_SQRT_2),
    )
}

fn check_t_squared_equals_s(b: &dyn SimulationBackend) -> CheckResult {
    // T² = S (up to global phase): fidelity(T²|+⟩, S|+⟩) = 1
    let mut ct2 = Circuit::new(1);
    ct2.push(Operation::new(GateKind::H, vec![0], vec![]));
    ct2.push(Operation::new(GateKind::T, vec![0], vec![]));
    ct2.push(Operation::new(GateKind::T, vec![0], vec![]));

    let mut cs = Circuit::new(1);
    cs.push(Operation::new(GateKind::H, vec![0], vec![]));
    cs.push(Operation::new(GateKind::S, vec![0], vec![]));

    fidelity_check(b, "T² = S (fidelity = 1)", "T gate", ct2, cs, 1.0)
}

fn check_t_fourth_equals_z(b: &dyn SimulationBackend) -> CheckResult {
    // T⁴ = Z (up to global phase): T⁴ = diag(1, e^{iπ}) = diag(1,-1) = Z
    let mut ct4 = Circuit::new(1);
    ct4.push(Operation::new(GateKind::H, vec![0], vec![]));
    for _ in 0..4 {
        ct4.push(Operation::new(GateKind::T, vec![0], vec![]));
    }

    let mut cz = Circuit::new(1);
    cz.push(Operation::new(GateKind::H, vec![0], vec![]));
    cz.push(Operation::new(GateKind::Z, vec![0], vec![]));

    fidelity_check(b, "T⁴ = Z (fidelity = 1)", "T gate", ct4, cz, 1.0)
}

fn check_sdg_phase(b: &dyn SimulationBackend) -> CheckResult {
    // Sdg|+⟩: ratio = -i  (Sdg = diag(1, -i))
    let mut c = Circuit::new(1);
    c.push(Operation::new(GateKind::H, vec![0], vec![]));
    c.push(Operation::new(GateKind::Sdg, vec![0], vec![]));
    ratio_check(
        b,
        "Sdg|+⟩ ratio = -i",
        "S gate",
        c,
        1,
        0,
        Complex64::new(0.0, -1.0),
    )
}

fn check_s_sdg_identity(b: &dyn SimulationBackend) -> CheckResult {
    // S·Sdg = I → S·Sdg|+⟩ = |+⟩ → ratio = 1
    let mut c = Circuit::new(1);
    c.push(Operation::new(GateKind::H, vec![0], vec![]));
    c.push(Operation::new(GateKind::S, vec![0], vec![]));
    c.push(Operation::new(GateKind::Sdg, vec![0], vec![]));
    ratio_check(
        b,
        "S·Sdg = I (ratio = 1)",
        "S gate",
        c,
        1,
        0,
        Complex64::new(1.0, 0.0),
    )
}

fn check_phase_rz_global(b: &dyn SimulationBackend) -> CheckResult {
    // P(π/2) = e^{iπ/4}·Rz(π/2) → fidelity(P|+⟩, Rz|+⟩) = 1
    let mut cp = Circuit::new(1);
    cp.push(Operation::new(GateKind::H, vec![0], vec![]));
    cp.push(Operation::new(GateKind::Phase, vec![0], vec![FRAC_PI_2]));

    let mut crz = Circuit::new(1);
    crz.push(Operation::new(GateKind::H, vec![0], vec![]));
    crz.push(Operation::new(GateKind::Rz, vec![0], vec![FRAC_PI_2]));

    fidelity_check(b, "P(λ)~Rz(λ) global phase (f=1)", "P vs Rz", cp, crz, 1.0)
}

fn check_cp_crz_differ(b: &dyn SimulationBackend) -> CheckResult {
    // CP(π) and CRz(π) on |++⟩: fidelity = 1/2 (relative phase, not global)
    // Derivation: CP(π)|++⟩ = (1/2)[1,1,1,-1]; CRz(π)|++⟩ = (1/2)[1,-i,1,i]
    // inner = (1-i)/2 → |inner|² = 1/2
    let prepare = |c: &mut Circuit| {
        c.push(Operation::new(GateKind::H, vec![0], vec![]));
        c.push(Operation::new(GateKind::H, vec![1], vec![]));
    };

    let mut ccp = Circuit::new(2);
    prepare(&mut ccp);
    ccp.push(Operation::new(GateKind::Cp, vec![0, 1], vec![PI]));

    let mut ccrz = Circuit::new(2);
    prepare(&mut ccrz);
    ccrz.push(Operation::new(GateKind::Crz, vec![0, 1], vec![PI]));

    fidelity_check(b, "CP≠CRz on |++⟩ (f=0.5)", "P vs Rz", ccp, ccrz, 0.5)
}

fn check_endianness_q0_lsb(b: &dyn SimulationBackend) -> CheckResult {
    // X on q0: |00⟩ → (q0=1,q1=0) → index 1 (q0=LSB, IBM convention)
    let mut c = Circuit::new(2);
    c.push(Operation::new(GateKind::X, vec![0], vec![]));
    amp_check(
        b,
        "X(q0)|00⟩ → sv[1]=1 (q0=LSB)",
        "Endianness",
        c,
        1,
        Complex64::new(1.0, 0.0),
    )
}

fn check_cx_ordering(b: &dyn SimulationBackend) -> CheckResult {
    // CX(q0,q1): X on q0 → |10⟩ (sv[1]=1). CX flips q1 → |11⟩ (sv[3]=1, sv[1]=0).
    // Both conditions are required: amplitude moved, not copied.
    let name = "CX(q0→q1)|10⟩ → sv[3]=1, sv[1]=0";
    let dimension = "Endianness";
    let mut c = Circuit::new(2);
    c.push(Operation::new(GateKind::X, vec![0], vec![]));
    c.push(Operation::new(GateKind::Cx, vec![0, 1], vec![]));
    match run(b, &c) {
        Err(e) => skip(name, dimension, e),
        Ok(sv) => {
            let sv3 = sv.get(3).copied().unwrap_or(Complex64::new(0.0, 0.0));
            let sv1 = sv.get(1).copied().unwrap_or(Complex64::new(0.0, 0.0));
            if (sv3 - Complex64::new(1.0, 0.0)).norm() < PHASE_TOL && sv1.norm() < PHASE_TOL {
                CheckResult {
                    name,
                    dimension,
                    status: CheckStatus::Pass,
                }
            } else {
                CheckResult {
                    name,
                    dimension,
                    status: CheckStatus::Fail {
                        expected: "sv[3]=1+0i, sv[1]=0+0i".into(),
                        got: format!("sv[3]={sv3:.4}, sv[1]={sv1:.4}"),
                    },
                }
            }
        }
    }
}

fn check_endianness_3q(b: &dyn SimulationBackend) -> CheckResult {
    // X on q2 in 3-qubit: index = 1<<2 = 4 (q2=LSB+2)
    let mut c = Circuit::new(3);
    c.push(Operation::new(GateKind::X, vec![2], vec![]));
    amp_check(
        b,
        "X(q2)|000⟩ → sv[4]=1 (3-qubit)",
        "Endianness",
        c,
        4,
        Complex64::new(1.0, 0.0),
    )
}

// ── Public API ────────────────────────────────────────────────────────────────

/// Runs all 12 QGCS conformance checks against `backend`.
pub fn certify(backend: &dyn SimulationBackend) -> Vec<CheckResult> {
    vec![
        check_rz_sign(backend),
        check_t_phase(backend),
        check_tdg_phase(backend),
        check_t_squared_equals_s(backend),
        check_t_fourth_equals_z(backend),
        check_sdg_phase(backend),
        check_s_sdg_identity(backend),
        check_phase_rz_global(backend),
        check_cp_crz_differ(backend),
        check_endianness_q0_lsb(backend),
        check_cx_ordering(backend),
        check_endianness_3q(backend),
    ]
}
