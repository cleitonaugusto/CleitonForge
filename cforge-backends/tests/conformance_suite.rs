//! Quantum Gate Convention Conformance Suite (QGCS v0.1)
//!
//! Phase-sensitive tests verifying NativeStateVectorBackend against the
//! OpenQASM 3 gate convention standard. Each test checks amplitude ratios,
//! not just probabilities — catching sign-convention bugs that are exactly
//! invisible to probability-only benchmarks such as Quantum Volume.
//!
//! ## Verification methodology
//!
//! For every gate G and witness state |ψ⟩, the expected amplitude ratio
//! sv[1]/sv[0] (or the exact amplitude sv[k]) is derived analytically from
//! the OpenQASM 3 matrix definition BEFORE the circuit is written. The
//! assert comes first; the circuit is evidence, not the spec.
//!
//! ## Qubit ordering (statevector index convention)
//!
//! Little-endian: qubit 0 is bit 0 (LSB).
//! index = Σ q_i · 2^i
//!
//! Two-qubit states (q1, q0):
//!   |00⟩ → index 0,  |01⟩ (q0=1) → index 1
//!   |10⟩ (q1=1) → index 2,  |11⟩ → index 3

use cforge_backends::{NativeStateVectorBackend, SimulationBackend};
use cforge_core::{Circuit, GateKind, Operation};
use num_complex::Complex64;
use std::f64::consts::{FRAC_1_SQRT_2, FRAC_PI_2, PI};

fn sv(c: &Circuit) -> Vec<Complex64> {
    NativeStateVectorBackend.run(c, 0, 0).unwrap().statevector
}

fn fidelity(a: &[Complex64], b: &[Complex64]) -> f64 {
    let inner: Complex64 = a.iter().zip(b).map(|(x, y)| x.conj() * y).sum();
    inner.norm_sqr()
}

// ── T gate ────────────────────────────────────────────────────────────────────

/// OpenQASM 3: T = diag(1, e^{iπ/4}).
///
/// Derivation: T|+⟩ = (T|0⟩ + T|1⟩)/√2 = (|0⟩ + e^{iπ/4}|1⟩)/√2
/// → sv[1]/sv[0] = e^{iπ/4} = cos(π/4) + i·sin(π/4) = 1/√2 + i/√2.
/// The wrong convention T = diag(1, e^{-iπ/4}) gives ratio 1/√2 - i/√2.
#[test]
fn t_gate_phase_openqasm3() {
    let expected = Complex64::new(FRAC_1_SQRT_2, FRAC_1_SQRT_2);

    let mut c = Circuit::new(1);
    c.push(Operation::new(GateKind::H, vec![0], vec![]));
    c.push(Operation::new(GateKind::T, vec![0], vec![]));

    let result = sv(&c);
    let ratio = result[1] / result[0];
    assert!(
        (ratio - expected).norm() < 1e-10,
        "T|+⟩ phase wrong: sv[1]/sv[0] = {ratio:.6}, expected e^{{iπ/4}} = {expected:.6}\n\
         OpenQASM 3: T = diag(1, e^{{iπ/4}})"
    );
}

/// OpenQASM 3: Tdg = T† = diag(1, e^{-iπ/4}).
///
/// Derivation: Tdg|+⟩: ratio = e^{-iπ/4} = 1/√2 - i/√2.
#[test]
fn tdg_gate_phase_openqasm3() {
    let expected = Complex64::new(FRAC_1_SQRT_2, -FRAC_1_SQRT_2);

    let mut c = Circuit::new(1);
    c.push(Operation::new(GateKind::H, vec![0], vec![]));
    c.push(Operation::new(GateKind::Tdg, vec![0], vec![]));

    let result = sv(&c);
    let ratio = result[1] / result[0];
    assert!(
        (ratio - expected).norm() < 1e-10,
        "Tdg|+⟩ phase wrong: ratio = {ratio:.6}, expected e^{{-iπ/4}} = {expected:.6}"
    );
}

/// T applied twice must equal S up to global phase.
///
/// Derivation: T² = diag(1, (e^{iπ/4})²) = diag(1, e^{iπ/2}) = diag(1, i) = S.
/// Fidelity(T²|+⟩, S|+⟩) = 1.
#[test]
fn t_squared_equals_s_gate() {
    let mut c_t2 = Circuit::new(1);
    c_t2.push(Operation::new(GateKind::H, vec![0], vec![]));
    c_t2.push(Operation::new(GateKind::T, vec![0], vec![]));
    c_t2.push(Operation::new(GateKind::T, vec![0], vec![]));

    let mut c_s = Circuit::new(1);
    c_s.push(Operation::new(GateKind::H, vec![0], vec![]));
    c_s.push(Operation::new(GateKind::S, vec![0], vec![]));

    let f = fidelity(&sv(&c_t2), &sv(&c_s));
    assert!(
        (f - 1.0).abs() < 1e-10,
        "T² must equal S up to global phase; fidelity = {f}"
    );
}

/// T and Tdg are mutual inverses: T†T = I.
///
/// Derivation: T·Tdg = diag(1, e^{iπ/4})·diag(1, e^{-iπ/4}) = I.
/// T†T|+⟩ = |+⟩ → ratio = 1.
#[test]
fn t_times_tdg_is_identity() {
    let mut c = Circuit::new(1);
    c.push(Operation::new(GateKind::H, vec![0], vec![]));
    c.push(Operation::new(GateKind::T, vec![0], vec![]));
    c.push(Operation::new(GateKind::Tdg, vec![0], vec![]));

    let result = sv(&c);
    let ratio = result[1] / result[0];
    assert!(
        (ratio - Complex64::new(1.0, 0.0)).norm() < 1e-10,
        "T·Tdg|+⟩ must return |+⟩; ratio = {ratio:.6}, expected 1"
    );
}

// ── S gate ────────────────────────────────────────────────────────────────────

/// OpenQASM 3: Sdg = S† = diag(1, -i) = diag(1, e^{-iπ/2}).
///
/// Derivation: Sdg|+⟩: ratio = e^{-iπ/2} = -i.
/// S = diag(1, i) already tested via rz_pi2_equivalent_to_s_gate in equivalence.rs.
#[test]
fn sdg_gate_phase_openqasm3() {
    let expected = Complex64::new(0.0, -1.0);

    let mut c = Circuit::new(1);
    c.push(Operation::new(GateKind::H, vec![0], vec![]));
    c.push(Operation::new(GateKind::Sdg, vec![0], vec![]));

    let result = sv(&c);
    let ratio = result[1] / result[0];
    assert!(
        (ratio - expected).norm() < 1e-10,
        "Sdg|+⟩ phase wrong: ratio = {ratio:.6}, expected -i"
    );
}

/// S · Sdg = I: applying S then Sdg recovers the original state.
///
/// Derivation: diag(1,i)·diag(1,-i) = diag(1,1) = I.
/// S·Sdg|+⟩ = |+⟩ → ratio = 1.
#[test]
fn s_times_sdg_is_identity() {
    let mut c = Circuit::new(1);
    c.push(Operation::new(GateKind::H, vec![0], vec![]));
    c.push(Operation::new(GateKind::S, vec![0], vec![]));
    c.push(Operation::new(GateKind::Sdg, vec![0], vec![]));

    let result = sv(&c);
    let ratio = result[1] / result[0];
    assert!(
        (ratio - Complex64::new(1.0, 0.0)).norm() < 1e-10,
        "S·Sdg|+⟩ must return |+⟩; ratio = {ratio:.6}, expected 1"
    );
}

/// T⁴ = S² = Z: four T gates must equal Z up to global phase.
///
/// Derivation: T⁴ = diag(1, e^{iπ}) = diag(1, -1) = Z.
/// Fidelity(T⁴|+⟩, Z|+⟩) = 1.
#[test]
fn t_fourth_equals_z_gate() {
    let mut c_t4 = Circuit::new(1);
    c_t4.push(Operation::new(GateKind::H, vec![0], vec![]));
    for _ in 0..4 {
        c_t4.push(Operation::new(GateKind::T, vec![0], vec![]));
    }

    let mut c_z = Circuit::new(1);
    c_z.push(Operation::new(GateKind::H, vec![0], vec![]));
    c_z.push(Operation::new(GateKind::Z, vec![0], vec![]));

    let f = fidelity(&sv(&c_t4), &sv(&c_z));
    assert!(
        (f - 1.0).abs() < 1e-10,
        "T⁴ must equal Z up to global phase; fidelity = {f}"
    );
}

// ── Phase gate P(λ) vs Rz(λ) ─────────────────────────────────────────────────

/// P(λ) = diag(1, e^{iλ}) and Rz(λ) = diag(e^{-iλ/2}, e^{+iλ/2}) differ
/// only by the global phase e^{iλ/2}: P(λ) = e^{iλ/2}·Rz(λ).
///
/// Derivation for λ = π/2:
///   P(π/2)|+⟩: sv = [1/√2, i/√2], ratio = i
///   Rz(π/2)|+⟩: sv = [e^{-iπ/4}/√2, e^{+iπ/4}/√2], ratio = e^{iπ/2} = i
///   Inner product = (1/2)(e^{-iπ/4} + i·e^{+iπ/4}·(−i)) ... = e^{iπ/4}
///   |inner|² = 1 → fidelity = 1.
#[test]
fn phase_gate_equivalent_to_rz_up_to_global_phase() {
    let lambda = FRAC_PI_2;

    let mut c_p = Circuit::new(1);
    c_p.push(Operation::new(GateKind::H, vec![0], vec![]));
    c_p.push(Operation::new(GateKind::Phase, vec![0], vec![lambda]));

    let mut c_rz = Circuit::new(1);
    c_rz.push(Operation::new(GateKind::H, vec![0], vec![]));
    c_rz.push(Operation::new(GateKind::Rz, vec![0], vec![lambda]));

    let f = fidelity(&sv(&c_p), &sv(&c_rz));
    assert!(
        (f - 1.0).abs() < 1e-10,
        "P(π/2) and Rz(π/2) must agree up to global phase; fidelity = {f}"
    );
}

/// CP(λ) and CRz(λ) are NOT equivalent even up to global phase: when used as
/// controlled gates, the global-phase difference becomes a relative phase that
/// is measurable on superposition inputs.
///
/// Setup: |++⟩ = H⊗H|00⟩, control=q0, target=q1.
///
/// Derivation:
///   |++⟩ = (1/2)[1, 1, 1, 1]  (indices 0..3 in little-endian)
///
///   CP(π)|++⟩: P(π)|1⟩ = e^{iπ}|1⟩ = -|1⟩
///     only index 3 (q0=1, q1=1) is affected: multiply by -1
///     result = (1/2)[1, 1, 1, -1]  → sv[3] = -1/2
///
///   CRz(π)|++⟩: Rz(π) = diag(e^{-iπ/2}, e^{+iπ/2}) = diag(-i, i)
///     index 1 (q0=1, q1=0): target qubit |0⟩ → multiply by e^{-iπ/2} = -i  → sv[1] = -i/2
///     index 3 (q0=1, q1=1): target qubit |1⟩ → multiply by e^{+iπ/2} = i   → sv[3] = i/2
///     result = (1/2)[1, -i, 1, i]
///
///   Fidelity = |⟨CP|CRz⟩|²:
///     inner = (1/4)(1 + 1·(-i) + 1 + (-1)·i) = (1/4)(2 - 2i) = (1-i)/2
///     |inner|² = |(1-i)/2|² = 2/4 = 1/2
#[test]
fn cp_differs_from_crz_on_superposition() {
    let prepare = |c: &mut Circuit| {
        c.push(Operation::new(GateKind::H, vec![0], vec![]));
        c.push(Operation::new(GateKind::H, vec![1], vec![]));
    };

    let mut c_cp = Circuit::new(2);
    prepare(&mut c_cp);
    c_cp.push(Operation::new(GateKind::Cp, vec![0, 1], vec![PI]));
    let sv_cp = sv(&c_cp);

    let mut c_crz = Circuit::new(2);
    prepare(&mut c_crz);
    c_crz.push(Operation::new(GateKind::Crz, vec![0, 1], vec![PI]));
    let sv_crz = sv(&c_crz);

    assert!(
        (sv_cp[3] - Complex64::new(-0.5, 0.0)).norm() < 1e-10,
        "CP(π)|++⟩: sv[3] expected -1/2, got {:.6}",
        sv_cp[3]
    );
    assert!(
        (sv_crz[1] - Complex64::new(0.0, -0.5)).norm() < 1e-10,
        "CRz(π)|++⟩: sv[1] expected -i/2, got {:.6}",
        sv_crz[1]
    );
    assert!(
        (sv_crz[3] - Complex64::new(0.0, 0.5)).norm() < 1e-10,
        "CRz(π)|++⟩: sv[3] expected i/2, got {:.6}",
        sv_crz[3]
    );

    let f = fidelity(&sv_cp, &sv_crz);
    assert!(
        (f - 0.5).abs() < 1e-10,
        "CP and CRz must differ on |++⟩; fidelity = {f}, expected 0.5"
    );
}

// ── Endianness / qubit ordering ───────────────────────────────────────────────

/// OpenQASM 3 / IBM convention: q[0] is the least-significant bit (LSB).
///
/// Derivation (little-endian, index = q1·2 + q0):
///   X on q0: |00⟩ → q0 flips → (q0=1, q1=0) → index = 0·2 + 1 = 1 → sv[1] = 1
///   X on q1: |00⟩ → q1 flips → (q0=0, q1=1) → index = 1·2 + 0 = 2 → sv[2] = 1
#[test]
fn qubit_ordering_q0_is_lsb() {
    let mut c0 = Circuit::new(2);
    c0.push(Operation::new(GateKind::X, vec![0], vec![]));
    let sv0 = sv(&c0);
    assert!(
        (sv0[1] - Complex64::new(1.0, 0.0)).norm() < 1e-10,
        "X on q0: expected sv[1]=1 (q0=LSB), got sv[1]={:.6}, sv[2]={:.6}",
        sv0[1],
        sv0[2]
    );

    let mut c1 = Circuit::new(2);
    c1.push(Operation::new(GateKind::X, vec![1], vec![]));
    let sv1 = sv(&c1);
    assert!(
        (sv1[2] - Complex64::new(1.0, 0.0)).norm() < 1e-10,
        "X on q1: expected sv[2]=1 (q1 is bit 1), got sv[1]={:.6}, sv[2]={:.6}",
        sv1[1],
        sv1[2]
    );
}

/// CX semantics: qubits[0] is control, qubits[1] is target.
///
/// Derivation:
///   X on q0: |00⟩ → index 1 (q0=1, q1=0).
///   CX(ctrl=q0, target=q1): ctrl=1 triggers X on q1 → (q0=1, q1=1) → index 3.
///   sv[3] = 1.
#[test]
fn cx_control_target_ordering() {
    let mut c = Circuit::new(2);
    c.push(Operation::new(GateKind::X, vec![0], vec![]));
    c.push(Operation::new(GateKind::Cx, vec![0, 1], vec![]));
    let result = sv(&c);

    assert!(
        (result[3] - Complex64::new(1.0, 0.0)).norm() < 1e-10,
        "CX(q0→q1)|10⟩: expected |11⟩ at sv[3]=1, got sv[3]={:.6}",
        result[3]
    );
    assert!(
        result[1].norm() < 1e-10,
        "CX(q0→q1)|10⟩: sv[1] must be 0 after flip, got {:.6}",
        result[1]
    );
}

/// Three-qubit endianness: X on q2 maps |000⟩ → index 4 (bit 2 = value 4).
///
/// Derivation: index = q2·4 + q1·2 + q0·1
///   X on q2: (q0=0, q1=0, q2=1) → index = 4 → sv[4] = 1.
#[test]
fn qubit_ordering_3q_q2_is_bit2() {
    let mut c = Circuit::new(3);
    c.push(Operation::new(GateKind::X, vec![2], vec![]));
    let result = sv(&c);

    assert!(
        (result[4] - Complex64::new(1.0, 0.0)).norm() < 1e-10,
        "X on q2 in 3-qubit system: expected sv[4]=1, got sv[4]={:.6}",
        result[4]
    );
}
