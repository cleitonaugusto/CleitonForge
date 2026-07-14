//! Integration tests: same circuit through both backends must give
//! equivalent results — fidelity ≥ 0.9999 for all standard circuits.

use cforge_backends::{NativeStateVectorBackend, QuantRS2Backend, SimulationBackend};
use cforge_core::{Circuit, GateKind, Operation};

fn fidelity(a: &[num_complex::Complex64], b: &[num_complex::Complex64]) -> f64 {
    let inner: num_complex::Complex64 = a.iter().zip(b).map(|(x, y)| x.conj() * y).sum();
    inner.norm_sqr()
}

fn run_both(circuit: &Circuit) -> (Vec<num_complex::Complex64>, Vec<num_complex::Complex64>) {
    let r1 = NativeStateVectorBackend.run(circuit, 0, 0).unwrap();
    let r2 = QuantRS2Backend.run(circuit, 0, 0).unwrap();
    (r1.statevector, r2.statevector)
}

#[test]
fn bell_state_backends_agree() {
    let mut c = Circuit::new(2);
    c.push(Operation::new(GateKind::H, vec![0], vec![]));
    c.push(Operation::new(GateKind::Cx, vec![0, 1], vec![]));
    let (sv1, sv2) = run_both(&c);
    assert!(
        fidelity(&sv1, &sv2) > 0.9999,
        "fidelity = {}",
        fidelity(&sv1, &sv2)
    );
}

#[test]
fn ghz_state_backends_agree() {
    // 3-qubit GHZ: H q[0]; CX q[0],q[1]; CX q[0],q[2]
    let mut c = Circuit::new(3);
    c.push(Operation::new(GateKind::H, vec![0], vec![]));
    c.push(Operation::new(GateKind::Cx, vec![0, 1], vec![]));
    c.push(Operation::new(GateKind::Cx, vec![0, 2], vec![]));
    let (sv1, sv2) = run_both(&c);
    assert!(fidelity(&sv1, &sv2) > 0.9999);
}

#[test]
fn toffoli_backends_agree() {
    // Toffoli with |11⟩ control: X q[0]; X q[1]; CCX q[0],q[1],q[2]
    let mut c = Circuit::new(3);
    c.push(Operation::new(GateKind::X, vec![0], vec![]));
    c.push(Operation::new(GateKind::X, vec![1], vec![]));
    c.push(Operation::new(GateKind::Ccx, vec![0, 1, 2], vec![]));
    let (sv1, sv2) = run_both(&c);
    assert!(fidelity(&sv1, &sv2) > 0.9999);
}

#[test]
fn pauli_gates_backends_agree() {
    let mut c = Circuit::new(2);
    c.push(Operation::new(GateKind::X, vec![0], vec![]));
    c.push(Operation::new(GateKind::Y, vec![1], vec![]));
    c.push(Operation::new(GateKind::Z, vec![0], vec![]));
    let (sv1, sv2) = run_both(&c);
    assert!(fidelity(&sv1, &sv2) > 0.9999);
}

#[test]
fn rx_gate_backends_agree() {
    let angle = std::f64::consts::FRAC_PI_4;
    let mut c = Circuit::new(1);
    c.push(Operation::new(GateKind::Rx, vec![0], vec![angle]));
    let (sv1, sv2) = run_both(&c);
    assert!(fidelity(&sv1, &sv2) > 0.9999);
}

/// Rz convention divergence: quantrs2 applies [[e^{+iθ/2},0],[0,e^{-iθ/2}]] while
/// native statevector applies [[e^{-iθ/2},0],[0,e^{+iθ/2}]]. On a |0⟩ or |1⟩
/// input the difference is a global phase (fidelity = 1). On a superposition
/// the two components accumulate opposite relative phases, making states
/// physically distinct (fidelity < 1). This is documented in quantrs2_backend.rs
/// and is precisely the kind of inter-framework divergence CleitonForge surfaces.
#[test]
fn rz_on_basis_state_is_global_phase_only() {
    // Rz(θ)|0⟩ differs only by global phase between backends → fidelity = 1.
    use std::f64::consts::PI;
    for &angle in &[PI / 4.0, -3.0 * PI / 2.0, PI, -PI / 3.0] {
        let mut c = Circuit::new(1);
        c.push(Operation::new(GateKind::Rz, vec![0], vec![angle]));
        let (sv1, sv2) = run_both(&c);
        let f = fidelity(&sv1, &sv2);
        assert!(f > 0.9999, "angle={angle}: fidelity = {f}");
    }
}

#[test]
fn rz_on_superposition_exposes_convention_divergence() {
    // H|0⟩ = |+⟩, then Rz(θ). The relative phase between |0⟩ and |1⟩ amplitudes
    // differs between backends — this is the documented quantrs2 Rz divergence.
    let angle = std::f64::consts::FRAC_PI_4;
    let mut c = Circuit::new(1);
    c.push(Operation::new(GateKind::H, vec![0], vec![]));
    c.push(Operation::new(GateKind::Rz, vec![0], vec![angle]));
    let (sv1, sv2) = run_both(&c);
    let f = fidelity(&sv1, &sv2);
    // Fidelity should be < 1 — backends diverge on superposition inputs.
    assert!(
        f < 0.9999,
        "expected fidelity < 1 for Rz on superposition, got {f}"
    );
}

#[test]
fn grover_3q_backends_agree() {
    // 3-qubit Grover, target |101⟩, 2 iterations. Uses only H/X/CCX — no Rz.
    let mut c = Circuit::new(3);
    for q in 0..3 {
        c.push(Operation::new(GateKind::H, vec![q], vec![]));
    }
    for _ in 0..2 {
        c.push(Operation::new(GateKind::X, vec![1], vec![]));
        c.push(Operation::new(GateKind::H, vec![2], vec![]));
        c.push(Operation::new(GateKind::Ccx, vec![0, 1, 2], vec![]));
        c.push(Operation::new(GateKind::H, vec![2], vec![]));
        c.push(Operation::new(GateKind::X, vec![1], vec![]));
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
    let (sv1, sv2) = run_both(&c);
    assert!(
        fidelity(&sv1, &sv2) > 0.9999,
        "fidelity = {}",
        fidelity(&sv1, &sv2)
    );
}

/// Rz(π/2)|+⟩ = (|0⟩ + i|1⟩)/√2 — OpenQASM 3 phase conformance test.
///
/// OpenQASM 3: Rz(λ) = diag(e^{-iλ/2}, e^{+iλ/2}) → Rz(π/2) ~ S = diag(1, i).
/// A sign-inverted Rz gives ratio -i instead of +i — undetectable by
/// probability-only benchmarks (both yield 50/50 outcomes; QV is exactly blind).
#[test]
fn rz_phase_conformance_native_openqasm3() {
    use num_complex::Complex64;
    use std::f64::consts::FRAC_PI_2;

    let mut c = Circuit::new(1);
    c.push(Operation::new(GateKind::H, vec![0], vec![]));
    c.push(Operation::new(GateKind::Rz, vec![0], vec![FRAC_PI_2]));

    let sv = NativeStateVectorBackend.run(&c, 0, 0).unwrap().statevector;

    // sv[1]/sv[0] = e^{+iπ/4} / e^{-iπ/4} = e^{iπ/2} = i
    let ratio = sv[1] / sv[0];
    let diff = (ratio - Complex64::new(0.0, 1.0)).norm();
    assert!(
        diff < 1e-10,
        "Rz(π/2)|+⟩ phase wrong: sv[1]/sv[0] = {ratio:.6}, expected i\n\
         OpenQASM 3 requires Rz(π/2) = diag(e^{{-iπ/4}}, e^{{+iπ/4}})"
    );
}

/// Rz(π/2) must produce the same state as S = diag(1, i) up to global phase.
/// Fidelity |⟨ψ_rz|ψ_s⟩|² = 1 encodes the S-equivalence without phase sensitivity.
#[test]
fn rz_pi2_equivalent_to_s_gate() {
    use std::f64::consts::FRAC_PI_2;

    let mut c_rz = Circuit::new(1);
    c_rz.push(Operation::new(GateKind::H, vec![0], vec![]));
    c_rz.push(Operation::new(GateKind::Rz, vec![0], vec![FRAC_PI_2]));
    let sv_rz = NativeStateVectorBackend.run(&c_rz, 0, 0).unwrap().statevector;

    let mut c_s = Circuit::new(1);
    c_s.push(Operation::new(GateKind::H, vec![0], vec![]));
    c_s.push(Operation::new(GateKind::S, vec![0], vec![]));
    let sv_s = NativeStateVectorBackend.run(&c_s, 0, 0).unwrap().statevector;

    let f = fidelity(&sv_rz, &sv_s);
    assert!(
        (f - 1.0).abs() < 1e-10,
        "Rz(π/2) and S must be equivalent up to global phase; fidelity = {f}"
    );
}

/// quantrs2 Rz is inverted: diag(e^{+iθ/2}, e^{-iθ/2}) → ratio = -i for Rz(π/2)|+⟩.
/// This documents the known upstream bug. If quantrs2 ever fixes their convention,
/// this test will fail — update quantrs2_backend.rs compensator accordingly.
#[test]
fn quantrs2_rz_inverted_sign_documented() {
    use num_complex::Complex64;
    use std::f64::consts::FRAC_PI_2;

    let mut c = Circuit::new(1);
    c.push(Operation::new(GateKind::H, vec![0], vec![]));
    c.push(Operation::new(GateKind::Rz, vec![0], vec![FRAC_PI_2]));

    let sv = QuantRS2Backend.run(&c, 0, 0).unwrap().statevector;

    // quantrs2: diag(e^{+iπ/4}, e^{-iπ/4}) → ratio = e^{-iπ/2} = -i
    let ratio = sv[1] / sv[0];
    let diff = (ratio - Complex64::new(0.0, -1.0)).norm();
    assert!(
        diff < 1e-10,
        "quantrs2 Rz(π/2)|+⟩: expected ratio -i (known inverted convention), got {ratio:.6}"
    );
}
