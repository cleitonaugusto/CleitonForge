//! Quantum simulation backend backed by the `quantrs2-core` crate.
//!
//! Uses QuantRS2's gate matrix definitions as the source of truth for
//! unitary matrices and applies them to a state vector using
//! CleitonForge's own application algorithm. This proves that the
//! `SimulationBackend` trait architecture supports multiple independent
//! backends on the same circuit IR.
//!
//! Gate matrix conventions in `quantrs2-core` may differ from other
//! frameworks for some parametric gates (e.g. Rz phase convention).
//! Such discrepancies, when surfaced by running the same circuit
//! through both backends and comparing metrics, are precisely the kind
//! of inter-framework insight CleitonForge is designed to expose.

use std::collections::HashMap;

use num_complex::Complex64;
use quantrs2_core::{
    gate::functions::{
        single::{
            Hadamard, Identity, PGate, PauliX, PauliY, PauliZ, Phase, PhaseDagger, RotationX,
            RotationY, RotationZ, SqrtX, SqrtXDagger, TDagger, UGate, T,
        },
        GateOp,
    },
    qubit::QubitId,
};

use cforge_core::{Circuit, GateKind};

use crate::sample::sample_counts;
use crate::trait_def::{BackendError, SimulationBackend, SimulationResult};

const MAX_QUBITS: usize = 22;

/// Simulation backend that sources gate matrices from `quantrs2-core`.
pub struct QuantRS2Backend;

impl SimulationBackend for QuantRS2Backend {
    fn name(&self) -> &str {
        "statevector-quantrs2"
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
                "quantrs2 backend supports up to {MAX_QUBITS} qubits, circuit has {n}"
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

// ── Dummy QubitId for struct construction (matrix does not depend on it) ──

const Q0: QubitId = QubitId(0);

fn apply_gate(
    sv: &mut [Complex64],
    kind: GateKind,
    qubits: &[usize],
    params: &[f64],
) -> Result<(), String> {
    match kind {
        GateKind::Id => {
            let mat = Identity { target: Q0 }
                .matrix()
                .map_err(|e| e.to_string())?;
            apply1(sv, qubits[0], &to_u2(&mat));
        }
        GateKind::X => {
            let mat = PauliX { target: Q0 }.matrix().map_err(|e| e.to_string())?;
            apply1(sv, qubits[0], &to_u2(&mat));
        }
        GateKind::Y => {
            let mat = PauliY { target: Q0 }.matrix().map_err(|e| e.to_string())?;
            apply1(sv, qubits[0], &to_u2(&mat));
        }
        GateKind::Z => {
            let mat = PauliZ { target: Q0 }.matrix().map_err(|e| e.to_string())?;
            apply1(sv, qubits[0], &to_u2(&mat));
        }
        GateKind::H => {
            let mat = Hadamard { target: Q0 }
                .matrix()
                .map_err(|e| e.to_string())?;
            apply1(sv, qubits[0], &to_u2(&mat));
        }
        GateKind::S => {
            let mat = Phase { target: Q0 }.matrix().map_err(|e| e.to_string())?;
            apply1(sv, qubits[0], &to_u2(&mat));
        }
        GateKind::Sdg => {
            let mat = PhaseDagger { target: Q0 }
                .matrix()
                .map_err(|e| e.to_string())?;
            apply1(sv, qubits[0], &to_u2(&mat));
        }
        GateKind::T => {
            let mat = T { target: Q0 }.matrix().map_err(|e| e.to_string())?;
            apply1(sv, qubits[0], &to_u2(&mat));
        }
        GateKind::Tdg => {
            let mat = TDagger { target: Q0 }.matrix().map_err(|e| e.to_string())?;
            apply1(sv, qubits[0], &to_u2(&mat));
        }
        GateKind::Sx => {
            let mat = SqrtX { target: Q0 }.matrix().map_err(|e| e.to_string())?;
            apply1(sv, qubits[0], &to_u2(&mat));
        }
        GateKind::Sxdg => {
            let mat = SqrtXDagger { target: Q0 }
                .matrix()
                .map_err(|e| e.to_string())?;
            apply1(sv, qubits[0], &to_u2(&mat));
        }
        GateKind::Rx => {
            let mat = RotationX {
                target: Q0,
                theta: params[0],
            }
            .matrix()
            .map_err(|e| e.to_string())?;
            apply1(sv, qubits[0], &to_u2(&mat));
        }
        GateKind::Ry => {
            let mat = RotationY {
                target: Q0,
                theta: params[0],
            }
            .matrix()
            .map_err(|e| e.to_string())?;
            apply1(sv, qubits[0], &to_u2(&mat));
        }
        GateKind::Rz => {
            let mat = RotationZ {
                target: Q0,
                theta: params[0],
            }
            .matrix()
            .map_err(|e| e.to_string())?;
            apply1(sv, qubits[0], &to_u2(&mat));
        }
        GateKind::Phase => {
            let mat = PGate {
                target: Q0,
                lambda: params[0],
            }
            .matrix()
            .map_err(|e| e.to_string())?;
            apply1(sv, qubits[0], &to_u2(&mat));
        }
        GateKind::U => {
            let mat = UGate {
                target: Q0,
                theta: params[0],
                phi: params[1],
                lambda: params[2],
            }
            .matrix()
            .map_err(|e| e.to_string())?;
            apply1(sv, qubits[0], &to_u2(&mat));
        }
        GateKind::Cx => {
            let mat = PauliX { target: Q0 }.matrix().map_err(|e| e.to_string())?;
            apply_controlled1(sv, qubits[0], qubits[1], &to_u2(&mat));
        }
        GateKind::Cy => {
            let mat = PauliY { target: Q0 }.matrix().map_err(|e| e.to_string())?;
            apply_controlled1(sv, qubits[0], qubits[1], &to_u2(&mat));
        }
        GateKind::Cz => {
            let mat = PauliZ { target: Q0 }.matrix().map_err(|e| e.to_string())?;
            apply_controlled1(sv, qubits[0], qubits[1], &to_u2(&mat));
        }
        GateKind::Ch => {
            let mat = Hadamard { target: Q0 }
                .matrix()
                .map_err(|e| e.to_string())?;
            apply_controlled1(sv, qubits[0], qubits[1], &to_u2(&mat));
        }
        GateKind::Csx => {
            let mat = SqrtX { target: Q0 }.matrix().map_err(|e| e.to_string())?;
            apply_controlled1(sv, qubits[0], qubits[1], &to_u2(&mat));
        }
        GateKind::Crx => {
            let mat = RotationX {
                target: Q0,
                theta: params[0],
            }
            .matrix()
            .map_err(|e| e.to_string())?;
            apply_controlled1(sv, qubits[0], qubits[1], &to_u2(&mat));
        }
        GateKind::Cry => {
            let mat = RotationY {
                target: Q0,
                theta: params[0],
            }
            .matrix()
            .map_err(|e| e.to_string())?;
            apply_controlled1(sv, qubits[0], qubits[1], &to_u2(&mat));
        }
        GateKind::Crz => {
            let mat = RotationZ {
                target: Q0,
                theta: params[0],
            }
            .matrix()
            .map_err(|e| e.to_string())?;
            apply_controlled1(sv, qubits[0], qubits[1], &to_u2(&mat));
        }
        GateKind::Cp => {
            let mat = PGate {
                target: Q0,
                lambda: params[0],
            }
            .matrix()
            .map_err(|e| e.to_string())?;
            apply_controlled1(sv, qubits[0], qubits[1], &to_u2(&mat));
        }
        GateKind::Cu => {
            let mat = UGate {
                target: Q0,
                theta: params[0],
                phi: params[1],
                lambda: params[2],
            }
            .matrix()
            .map_err(|e| e.to_string())?;
            apply_controlled1(sv, qubits[0], qubits[1], &to_u2(&mat));
        }
        GateKind::Swap => apply_swap(sv, qubits[0], qubits[1]),
        GateKind::Ccx => apply_ccx(sv, qubits[0], qubits[1], qubits[2]),
        GateKind::Cswap => apply_cswap(sv, qubits[0], qubits[1], qubits[2]),
    }
    Ok(())
}

// ── Matrix format conversion ──────────────────────────────────────────────────

type U2 = [[Complex64; 2]; 2];

fn to_u2(flat: &[Complex64]) -> U2 {
    [[flat[0], flat[1]], [flat[2], flat[3]]]
}

// ── State-vector application (same algorithm as native backend) ───────────────

fn apply1(sv: &mut [Complex64], k: usize, u: &U2) {
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

fn apply_controlled1(sv: &mut [Complex64], ctrl: usize, target: usize, u: &U2) {
    let ctrl_bit = 1 << ctrl;
    let tgt_bit = 1 << target;
    for i in 0..sv.len() {
        if (i & ctrl_bit != 0) && (i & tgt_bit == 0) {
            let i1 = i | tgt_bit;
            let a = sv[i];
            let b = sv[i1];
            sv[i] = u[0][0] * a + u[0][1] * b;
            sv[i1] = u[1][0] * a + u[1][1] * b;
        }
    }
}

fn apply_swap(sv: &mut [Complex64], a: usize, b: usize) {
    let bit_a = 1 << a;
    let bit_b = 1 << b;
    for i in 0..sv.len() {
        if (i & bit_a == 0) && (i & bit_b != 0) {
            sv.swap(i, (i | bit_a) & !bit_b);
        }
    }
}

fn apply_ccx(sv: &mut [Complex64], ctrl0: usize, ctrl1: usize, target: usize) {
    let cb0 = 1 << ctrl0;
    let cb1 = 1 << ctrl1;
    let tgt = 1 << target;
    for i in 0..sv.len() {
        if (i & cb0 != 0) && (i & cb1 != 0) && (i & tgt == 0) {
            sv.swap(i, i | tgt);
        }
    }
}

fn apply_cswap(sv: &mut [Complex64], ctrl: usize, a: usize, b: usize) {
    let cb = 1 << ctrl;
    let bit_a = 1 << a;
    let bit_b = 1 << b;
    for i in 0..sv.len() {
        if (i & cb != 0) && (i & bit_a == 0) && (i & bit_b != 0) {
            sv.swap(i, (i | bit_a) & !bit_b);
        }
    }
}
