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
    assert!(fidelity(&sv1, &sv2) > 0.9999, "fidelity = {}", fidelity(&sv1, &sv2));
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
