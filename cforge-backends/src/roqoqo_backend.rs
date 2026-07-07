//! Simulation backend sourcing gate matrices from the `roqoqo` crate (HQS).
//!
//! Uses roqoqo's unitary matrix definitions as the source of truth and applies
//! them to a statevector using CleitonForge's own application algorithm.
//!
//! **Rz convention**: roqoqo uses [[e^{-iθ/2},0],[0,e^{+iθ/2}]], which matches
//! CleitonForge's native statevector backend (IBM/Qiskit standard). This
//! differs from quantrs2-core, which reverses the sign. Running the same QAOA
//! circuit through all three backends surfaces this divergence explicitly —
//! roqoqo and native agree; quantrs2 gives the opposite cut states.
//!
//! **Matrix ordering**: 2-qubit matrices use the basis {|00⟩,|01⟩,|10⟩,|11⟩}
//! where q0 (first qubit argument) is the more-significant bit in the
//! subspace index. This matches our `apply2` implementation below.

use std::collections::HashMap;

use ndarray::Array2;
use num_complex::Complex64;
use roqoqo::operations::{
    ControlledPauliY, ControlledPauliZ, ControlledPhaseShift, Hadamard, InvSGate, InvSqrtPauliX,
    PauliX, PauliY, PauliZ, PhaseShiftState1, RotateX, RotateY, RotateZ, SGate, SqrtPauliX, TGate,
    Toffoli, CNOT, SWAP,
};
use roqoqo::prelude::OperateGate;

use cforge_core::{Circuit, GateKind};

use crate::sample::sample_counts;
use crate::trait_def::{BackendError, SimulationBackend, SimulationResult};

const MAX_QUBITS: usize = 22;

/// Simulation backend that sources gate matrices from `roqoqo` (HQS).
pub struct RoqoqoBackend;

impl SimulationBackend for RoqoqoBackend {
    fn name(&self) -> &str {
        "statevector-roqoqo"
    }

    fn run(
        &self,
        circuit: &Circuit,
        shots: usize,
        seed: u64,
    ) -> Result<SimulationResult, BackendError> {
        let n = circuit.num_qubits();
        if n > MAX_QUBITS {
            return Err(BackendError(format!(
                "roqoqo backend supports up to {MAX_QUBITS} qubits, circuit has {n}"
            )));
        }

        let mut sv = vec![Complex64::new(0.0, 0.0); 1 << n];
        sv[0] = Complex64::new(1.0, 0.0);

        for (op_idx, op) in circuit.operations.iter().enumerate() {
            apply_gate(&mut sv, op.kind, &op.qubits, &op.params)
                .map_err(|e| BackendError(format!("operation {op_idx}: {e}")))?;
        }

        let counts = if shots > 0 {
            sample_counts(&sv, shots, seed)
        } else {
            HashMap::new()
        };

        Ok(SimulationResult {
            statevector: sv,
            counts,
        })
    }
}

// ── Dummy qubit index for matrix-only construction ────────────────────────────

const Q0: usize = 0;
const Q1: usize = 1;
const Q2: usize = 2;

fn cf(v: f64) -> qoqo_calculator::CalculatorFloat {
    qoqo_calculator::CalculatorFloat::from(v)
}

// ── Gate dispatch ─────────────────────────────────────────────────────────────

fn apply_gate(
    sv: &mut [Complex64],
    kind: GateKind,
    qubits: &[usize],
    params: &[f64],
) -> Result<(), String> {
    let q = qubits;

    match kind {
        // ── Single-qubit non-parametric ──────────────────────────────────────
        GateKind::Id => {}
        GateKind::X => apply1(sv, q[0], mat1(PauliX::new(Q0).unitary_matrix())?),
        GateKind::Y => apply1(sv, q[0], mat1(PauliY::new(Q0).unitary_matrix())?),
        GateKind::Z => apply1(sv, q[0], mat1(PauliZ::new(Q0).unitary_matrix())?),
        GateKind::H => apply1(sv, q[0], mat1(Hadamard::new(Q0).unitary_matrix())?),
        GateKind::S => apply1(sv, q[0], mat1(SGate::new(Q0).unitary_matrix())?),
        GateKind::Sdg => apply1(sv, q[0], mat1(InvSGate::new(Q0).unitary_matrix())?),
        GateKind::T => apply1(sv, q[0], mat1(TGate::new(Q0).unitary_matrix())?),
        GateKind::Tdg => apply1(
            sv,
            q[0],
            conj_transpose_2x2(mat1(TGate::new(Q0).unitary_matrix())?),
        ),
        GateKind::Sx => apply1(sv, q[0], mat1(SqrtPauliX::new(Q0).unitary_matrix())?),
        GateKind::Sxdg => apply1(sv, q[0], mat1(InvSqrtPauliX::new(Q0).unitary_matrix())?),

        // ── Single-qubit parametric ──────────────────────────────────────────
        GateKind::Rx => apply1(
            sv,
            q[0],
            mat1(RotateX::new(Q0, cf(params[0])).unitary_matrix())?,
        ),
        GateKind::Ry => apply1(
            sv,
            q[0],
            mat1(RotateY::new(Q0, cf(params[0])).unitary_matrix())?,
        ),
        GateKind::Rz => apply1(
            sv,
            q[0],
            mat1(RotateZ::new(Q0, cf(params[0])).unitary_matrix())?,
        ),
        GateKind::Phase => apply1(
            sv,
            q[0],
            mat1(PhaseShiftState1::new(Q0, cf(params[0])).unitary_matrix())?,
        ),

        // U(θ, φ, λ): general single-qubit gate — computed directly.
        // U = [[cos(θ/2),          -e^{iλ}·sin(θ/2)],
        //      [e^{iφ}·sin(θ/2),   e^{i(φ+λ)}·cos(θ/2)]]
        GateKind::U => {
            let (th, phi, lam) = (params[0] / 2.0, params[1], params[2]);
            let c = th.cos();
            let s = th.sin();
            let ep = Complex64::from_polar(1.0, phi);
            let el = Complex64::from_polar(1.0, lam);
            let epl = Complex64::from_polar(1.0, phi + lam);
            apply1(
                sv,
                q[0],
                [[Complex64::new(c, 0.0), -el * s], [ep * s, epl * c]],
            );
        }

        // ── Two-qubit gates via roqoqo matrices ──────────────────────────────
        GateKind::Cx => apply2(sv, q[0], q[1], mat2(CNOT::new(Q0, Q1).unitary_matrix())?),
        GateKind::Cz => apply2(
            sv,
            q[0],
            q[1],
            mat2(ControlledPauliZ::new(Q0, Q1).unitary_matrix())?,
        ),
        GateKind::Cy => apply2(
            sv,
            q[0],
            q[1],
            mat2(ControlledPauliY::new(Q0, Q1).unitary_matrix())?,
        ),
        GateKind::Swap => apply2(sv, q[0], q[1], mat2(SWAP::new(Q0, Q1).unitary_matrix())?),
        GateKind::Cp => apply2(
            sv,
            q[0],
            q[1],
            mat2(ControlledPhaseShift::new(Q0, Q1, cf(params[0])).unitary_matrix())?,
        ),

        // Two-qubit gates without a direct roqoqo equivalent — use native math.
        GateKind::Ch
        | GateKind::Csx
        | GateKind::Crx
        | GateKind::Cry
        | GateKind::Crz
        | GateKind::Cu => apply_controlled_native(sv, kind, q, params),

        // ── Three-qubit gates ─────────────────────────────────────────────────
        GateKind::Ccx => {
            // Toffoli from roqoqo — apply as general 8×8 matrix.
            let u8 = mat3(Toffoli::new(Q0, Q1, Q2).unitary_matrix())?;
            apply3(sv, q[0], q[1], q[2], u8);
        }
        GateKind::Cswap => apply_cswap_native(sv, q[0], q[1], q[2]),
    }
    Ok(())
}

// ── Matrix format helpers ─────────────────────────────────────────────────────

type U2 = [[Complex64; 2]; 2];
type U4 = [[Complex64; 4]; 4];
type U8 = [[Complex64; 8]; 8];

fn mat1(r: Result<Array2<Complex64>, roqoqo::RoqoqoError>) -> Result<U2, String> {
    let m = r.map_err(|e| e.to_string())?;
    Ok([[m[[0, 0]], m[[0, 1]]], [m[[1, 0]], m[[1, 1]]]])
}

fn mat2(r: Result<Array2<Complex64>, roqoqo::RoqoqoError>) -> Result<U4, String> {
    let m = r.map_err(|e| e.to_string())?;
    let mut out = [[Complex64::new(0.0, 0.0); 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            out[i][j] = m[[i, j]];
        }
    }
    Ok(out)
}

fn mat3(r: Result<Array2<Complex64>, roqoqo::RoqoqoError>) -> Result<U8, String> {
    let m = r.map_err(|e| e.to_string())?;
    let mut out = [[Complex64::new(0.0, 0.0); 8]; 8];
    for i in 0..8 {
        for j in 0..8 {
            out[i][j] = m[[i, j]];
        }
    }
    Ok(out)
}

fn conj_transpose_2x2(u: U2) -> U2 {
    [
        [u[0][0].conj(), u[1][0].conj()],
        [u[0][1].conj(), u[1][1].conj()],
    ]
}

// ── Statevector application ───────────────────────────────────────────────────

fn apply1(sv: &mut [Complex64], k: usize, u: U2) {
    let stride = 1 << k;
    let len = sv.len();
    let mut i = 0;
    while i < len {
        for j in i..(i + stride) {
            let a = sv[j];
            let b = sv[j + stride];
            sv[j] = u[0][0] * a + u[0][1] * b;
            sv[j + stride] = u[1][0] * a + u[1][1] * b;
        }
        i += 2 * stride;
    }
}

/// Applies a 4×4 unitary to qubits (q0, q1).
///
/// Basis ordering: q0 is the more-significant bit in the 2-qubit subspace
/// index — matching roqoqo's CNOT(control=q0, target=q1) convention.
fn apply2(sv: &mut [Complex64], q0: usize, q1: usize, u: U4) {
    let bit0 = 1 << q0;
    let bit1 = 1 << q1;
    let len = sv.len();
    for i in 0..len {
        if (i & bit0 != 0) || (i & bit1 != 0) {
            continue;
        }
        // Four indices: subspace 0→i, 1→i|bit1, 2→i|bit0, 3→i|bit0|bit1
        let idx = [i, i | bit1, i | bit0, i | bit0 | bit1];
        let old: [Complex64; 4] = [sv[idx[0]], sv[idx[1]], sv[idx[2]], sv[idx[3]]];
        for row in 0..4 {
            sv[idx[row]] = (0..4).map(|col| u[row][col] * old[col]).sum();
        }
    }
}

/// Applies an 8×8 unitary to qubits (q0, q1, q2).
///
/// Subspace index: q0 is bit 2, q1 is bit 1, q2 is bit 0 of the 3-bit index.
fn apply3(sv: &mut [Complex64], q0: usize, q1: usize, q2: usize, u: U8) {
    let bit0 = 1 << q0;
    let bit1 = 1 << q1;
    let bit2 = 1 << q2;
    let len = sv.len();
    for i in 0..len {
        if (i & bit0 != 0) || (i & bit1 != 0) || (i & bit2 != 0) {
            continue;
        }
        let idx = [
            i,
            i | bit2,
            i | bit1,
            i | bit1 | bit2,
            i | bit0,
            i | bit0 | bit2,
            i | bit0 | bit1,
            i | bit0 | bit1 | bit2,
        ];
        let old: [Complex64; 8] = std::array::from_fn(|k| sv[idx[k]]);
        for row in 0..8 {
            sv[idx[row]] = (0..8).map(|col| u[row][col] * old[col]).sum();
        }
    }
}

// ── Native fallback for gates without roqoqo equivalents ─────────────────────

fn apply_controlled_native(sv: &mut [Complex64], kind: GateKind, q: &[usize], params: &[f64]) {
    let ctrl = q[0];
    let tgt = q[1];
    let cb = 1 << ctrl;
    let tb = 1 << tgt;

    let u: U2 = match kind {
        GateKind::Ch => {
            let s = 1.0 / std::f64::consts::SQRT_2;
            [
                [Complex64::new(s, 0.0), Complex64::new(s, 0.0)],
                [Complex64::new(s, 0.0), Complex64::new(-s, 0.0)],
            ]
        }
        GateKind::Csx => [
            [Complex64::new(0.5, 0.5), Complex64::new(0.5, -0.5)],
            [Complex64::new(0.5, -0.5), Complex64::new(0.5, 0.5)],
        ],
        GateKind::Crx => {
            let c = (params[0] / 2.0).cos();
            let s = (params[0] / 2.0).sin();
            [
                [Complex64::new(c, 0.0), Complex64::new(0.0, -s)],
                [Complex64::new(0.0, -s), Complex64::new(c, 0.0)],
            ]
        }
        GateKind::Cry => {
            let c = (params[0] / 2.0).cos();
            let s = (params[0] / 2.0).sin();
            [
                [Complex64::new(c, 0.0), Complex64::new(-s, 0.0)],
                [Complex64::new(s, 0.0), Complex64::new(c, 0.0)],
            ]
        }
        GateKind::Crz => {
            let h = params[0] / 2.0;
            [
                [Complex64::new(h.cos(), -h.sin()), Complex64::new(0.0, 0.0)],
                [Complex64::new(0.0, 0.0), Complex64::new(h.cos(), h.sin())],
            ]
        }
        GateKind::Cu => {
            let (th, phi, lam) = (params[0] / 2.0, params[1], params[2]);
            let ep = Complex64::from_polar(1.0, phi);
            let el = Complex64::from_polar(1.0, lam);
            let epl = Complex64::from_polar(1.0, phi + lam);
            [
                [Complex64::new(th.cos(), 0.0), -el * th.sin()],
                [ep * th.sin(), epl * th.cos()],
            ]
        }
        _ => return,
    };

    for i in 0..sv.len() {
        if (i & cb != 0) && (i & tb == 0) {
            let i1 = i | tb;
            let a = sv[i];
            let b = sv[i1];
            sv[i] = u[0][0] * a + u[0][1] * b;
            sv[i1] = u[1][0] * a + u[1][1] * b;
        }
    }
}

fn apply_cswap_native(sv: &mut [Complex64], ctrl: usize, a: usize, b: usize) {
    let cb = 1 << ctrl;
    let ba = 1 << a;
    let bb = 1 << b;
    for i in 0..sv.len() {
        if (i & cb != 0) && (i & ba == 0) && (i & bb != 0) {
            sv.swap(i, (i | ba) & !bb);
        }
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::trait_def::DEFAULT_SEED;
    use crate::{NativeStateVectorBackend, SimulationBackend};
    use cforge_core::Operation;

    fn fidelity(a: &[Complex64], b: &[Complex64]) -> f64 {
        let inner: Complex64 = a.iter().zip(b).map(|(x, y)| x.conj() * y).sum();
        inner.norm_sqr()
    }

    fn run_both(circuit: &Circuit) -> (Vec<Complex64>, Vec<Complex64>) {
        let r1 = NativeStateVectorBackend.run(circuit, 0, 0).unwrap();
        let r2 = RoqoqoBackend.run(circuit, 0, 0).unwrap();
        (r1.statevector, r2.statevector)
    }

    #[test]
    fn bell_state_roqoqo_agrees() {
        let mut c = Circuit::new(2);
        c.push(Operation::new(GateKind::H, vec![0], vec![]));
        c.push(Operation::new(GateKind::Cx, vec![0, 1], vec![]));
        let (sv1, sv2) = run_both(&c);
        let f = fidelity(&sv1, &sv2);
        assert!(f > 0.9999, "fidelity = {f}");
    }

    #[test]
    fn grover_3q_roqoqo_agrees() {
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
        let f = fidelity(&sv1, &sv2);
        assert!(f > 0.9999, "fidelity = {f}");
    }

    #[test]
    fn rz_on_superposition_roqoqo_agrees_with_native() {
        // roqoqo uses the same Rz convention as our native backend:
        // [[e^{-iθ/2}, 0], [0, e^{+iθ/2}]] — matching IBM/Qiskit standard.
        // This is the key difference from quantrs2 which reverses the sign.
        let angle = std::f64::consts::FRAC_PI_4;
        let mut c = Circuit::new(1);
        c.push(Operation::new(GateKind::H, vec![0], vec![]));
        c.push(Operation::new(GateKind::Rz, vec![0], vec![angle]));
        let (sv1, sv2) = run_both(&c);
        let f = fidelity(&sv1, &sv2);
        assert!(
            f > 0.9999,
            "roqoqo Rz diverges from native — expected agreement, got fidelity={f}"
        );
    }

    #[test]
    fn qaoa_maxcut_roqoqo_agrees_with_native() {
        // QAOA uses Rz → roqoqo should agree with statevector (both IBM convention).
        // quantrs2 gives the opposite cut states; roqoqo gives the correct ones.
        use std::f64::consts::PI;
        let gamma = -3.0 * PI / 4.0;
        let beta = -PI / 8.0;
        let mut c = Circuit::new(2);
        c.push(Operation::new(GateKind::H, vec![0], vec![]));
        c.push(Operation::new(GateKind::H, vec![1], vec![]));
        c.push(Operation::new(GateKind::Cx, vec![0, 1], vec![]));
        c.push(Operation::new(GateKind::Rz, vec![1], vec![2.0 * gamma]));
        c.push(Operation::new(GateKind::Cx, vec![0, 1], vec![]));
        c.push(Operation::new(GateKind::Rx, vec![0], vec![2.0 * beta]));
        c.push(Operation::new(GateKind::Rx, vec![1], vec![2.0 * beta]));
        let (sv1, sv2) = run_both(&c);
        let f = fidelity(&sv1, &sv2);
        assert!(f > 0.9999, "QAOA: roqoqo vs native fidelity = {f}");
    }

    #[test]
    fn roqoqo_backend_name() {
        assert_eq!(RoqoqoBackend.name(), "statevector-roqoqo");
    }

    #[test]
    fn roqoqo_shots_and_seed_deterministic() {
        let mut c = Circuit::new(2);
        c.push(Operation::new(GateKind::H, vec![0], vec![]));
        c.push(Operation::new(GateKind::Cx, vec![0, 1], vec![]));
        let r1 = RoqoqoBackend.run(&c, 1024, DEFAULT_SEED).unwrap();
        let r2 = RoqoqoBackend.run(&c, 1024, DEFAULT_SEED).unwrap();
        assert_eq!(r1.counts, r2.counts);
    }
}
