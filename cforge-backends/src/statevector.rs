//! Native state-vector backend written from scratch.
//!
//! Represents the quantum state as a `Vec<Complex64>` of length 2^n and
//! applies each gate as a unitary matrix multiplication on the
//! appropriate subset of amplitudes. Supports up to ~20 qubits (2^20
//! ≈ 10^6 complex entries, ~16 MB), sufficient for algorithm validation.

use std::collections::HashMap;
use std::f64::consts::{FRAC_1_SQRT_2, PI};

use num_complex::Complex64;

use cforge_core::{Circuit, GateKind};

use crate::sample::sample_counts;
use crate::trait_def::{BackendError, SimulationBackend, SimulationResult};

const MAX_QUBITS: usize = 22;

/// State-vector backend built entirely within CleitonForge.
pub struct NativeStateVectorBackend;

impl SimulationBackend for NativeStateVectorBackend {
    fn name(&self) -> &str {
        "statevector-native"
    }

    fn run(&self, circuit: &Circuit, shots: usize, seed: u64) -> Result<SimulationResult, BackendError> {
        let n = circuit.num_qubits();
        if n > MAX_QUBITS {
            return Err(BackendError(format!(
                "native backend supports up to {MAX_QUBITS} qubits, circuit has {n}"
            )));
        }

        let mut sv = vec![Complex64::new(0.0, 0.0); 1 << n];
        sv[0] = Complex64::new(1.0, 0.0); // |0...0⟩

        for (op_idx, op) in circuit.operations.iter().enumerate() {
            apply_gate(&mut sv, op.kind, &op.qubits, &op.params, n)
                .map_err(|e| BackendError(format!("operation {op_idx}: {e}")))?;
        }

        let counts = if shots > 0 {
            sample_counts(&sv, shots, seed)
        } else {
            HashMap::new()
        };

        Ok(SimulationResult { statevector: sv, counts })
    }
}

// ── Gate application ─────────────────────────────────────────────────────────

fn apply_gate(
    sv: &mut [Complex64],
    kind: GateKind,
    qubits: &[usize],
    params: &[f64],
    _n: usize,
) -> Result<(), String> {
    match kind {
        // Single-qubit gates with no parameters
        GateKind::Id => {}
        GateKind::X => apply1(sv, qubits[0], &gate_x()),
        GateKind::Y => apply1(sv, qubits[0], &gate_y()),
        GateKind::Z => apply1(sv, qubits[0], &gate_z()),
        GateKind::H => apply1(sv, qubits[0], &gate_h()),
        GateKind::S => apply1(sv, qubits[0], &gate_s()),
        GateKind::Sdg => apply1(sv, qubits[0], &gate_sdg()),
        GateKind::T => apply1(sv, qubits[0], &gate_t()),
        GateKind::Tdg => apply1(sv, qubits[0], &gate_tdg()),
        GateKind::Sx => apply1(sv, qubits[0], &gate_sx()),
        GateKind::Sxdg => apply1(sv, qubits[0], &gate_sxdg()),

        // Parametric single-qubit gates
        GateKind::Rx => apply1(sv, qubits[0], &gate_rx(params[0])),
        GateKind::Ry => apply1(sv, qubits[0], &gate_ry(params[0])),
        GateKind::Rz => apply1(sv, qubits[0], &gate_rz(params[0])),
        GateKind::Phase => apply1(sv, qubits[0], &gate_phase(params[0])),
        GateKind::U => apply1(sv, qubits[0], &gate_u(params[0], params[1], params[2])),

        // Two-qubit gates
        GateKind::Cx => apply_controlled1(sv, qubits[0], qubits[1], &gate_x()),
        GateKind::Cy => apply_controlled1(sv, qubits[0], qubits[1], &gate_y()),
        GateKind::Cz => apply_controlled1(sv, qubits[0], qubits[1], &gate_z()),
        GateKind::Ch => apply_controlled1(sv, qubits[0], qubits[1], &gate_h()),
        GateKind::Csx => apply_controlled1(sv, qubits[0], qubits[1], &gate_sx()),
        GateKind::Crx => apply_controlled1(sv, qubits[0], qubits[1], &gate_rx(params[0])),
        GateKind::Cry => apply_controlled1(sv, qubits[0], qubits[1], &gate_ry(params[0])),
        GateKind::Crz => apply_controlled1(sv, qubits[0], qubits[1], &gate_rz(params[0])),
        GateKind::Cp => apply_controlled1(sv, qubits[0], qubits[1], &gate_phase(params[0])),
        GateKind::Cu => apply_controlled1(
            sv,
            qubits[0],
            qubits[1],
            &gate_cu(params[0], params[1], params[2], params[3]),
        ),
        GateKind::Swap => apply_swap(sv, qubits[0], qubits[1]),

        // Three-qubit gates
        GateKind::Ccx => apply_ccx(sv, qubits[0], qubits[1], qubits[2]),
        GateKind::Cswap => apply_cswap(sv, qubits[0], qubits[1], qubits[2]),
    }
    Ok(())
}

// ── Single-qubit unitary application ────────────────────────────────────────

type U2 = [[Complex64; 2]; 2];

/// Applies a 2×2 unitary to qubit `k`.
fn apply1(sv: &mut [Complex64], k: usize, u: &U2) {
    let stride = 1 << k;
    let len = sv.len();
    let mut i = 0;
    while i < len {
        for j in i..(i + stride) {
            let i0 = j;
            let i1 = j + stride;
            let a = sv[i0];
            let b = sv[i1];
            sv[i0] = u[0][0] * a + u[0][1] * b;
            sv[i1] = u[1][0] * a + u[1][1] * b;
        }
        i += 2 * stride;
    }
}

/// Applies a 2×2 unitary to `target` qubit when `ctrl` qubit is |1⟩.
fn apply_controlled1(sv: &mut [Complex64], ctrl: usize, target: usize, u: &U2) {
    let ctrl_bit = 1 << ctrl;
    let tgt_bit = 1 << target;
    let len = sv.len();
    for i in 0..len {
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
    let len = sv.len();
    for i in 0..len {
        if (i & bit_a == 0) && (i & bit_b != 0) {
            let j = (i | bit_a) & !(bit_b);
            sv.swap(i, j);
        }
    }
}

fn apply_ccx(sv: &mut [Complex64], ctrl0: usize, ctrl1: usize, target: usize) {
    let cb0 = 1 << ctrl0;
    let cb1 = 1 << ctrl1;
    let tgt = 1 << target;
    let len = sv.len();
    for i in 0..len {
        if (i & cb0 != 0) && (i & cb1 != 0) && (i & tgt == 0) {
            sv.swap(i, i | tgt);
        }
    }
}

fn apply_cswap(sv: &mut [Complex64], ctrl: usize, a: usize, b: usize) {
    let cb = 1 << ctrl;
    let bit_a = 1 << a;
    let bit_b = 1 << b;
    let len = sv.len();
    for i in 0..len {
        if (i & cb != 0) && (i & bit_a == 0) && (i & bit_b != 0) {
            let j = (i | bit_a) & !bit_b;
            sv.swap(i, j);
        }
    }
}

// ── Gate matrix definitions ──────────────────────────────────────────────────

#[inline]
fn c(re: f64, im: f64) -> Complex64 {
    Complex64::new(re, im)
}

fn gate_x() -> U2 {
    [[c(0.0, 0.0), c(1.0, 0.0)], [c(1.0, 0.0), c(0.0, 0.0)]]
}
fn gate_y() -> U2 {
    [[c(0.0, 0.0), c(0.0, -1.0)], [c(0.0, 1.0), c(0.0, 0.0)]]
}
fn gate_z() -> U2 {
    [[c(1.0, 0.0), c(0.0, 0.0)], [c(0.0, 0.0), c(-1.0, 0.0)]]
}
fn gate_h() -> U2 {
    let f = FRAC_1_SQRT_2;
    [[c(f, 0.0), c(f, 0.0)], [c(f, 0.0), c(-f, 0.0)]]
}
fn gate_s() -> U2 {
    [[c(1.0, 0.0), c(0.0, 0.0)], [c(0.0, 0.0), c(0.0, 1.0)]]
}
fn gate_sdg() -> U2 {
    [[c(1.0, 0.0), c(0.0, 0.0)], [c(0.0, 0.0), c(0.0, -1.0)]]
}
fn gate_t() -> U2 {
    let e = Complex64::from_polar(1.0, PI / 4.0);
    [[c(1.0, 0.0), c(0.0, 0.0)], [c(0.0, 0.0), e]]
}
fn gate_tdg() -> U2 {
    let e = Complex64::from_polar(1.0, -PI / 4.0);
    [[c(1.0, 0.0), c(0.0, 0.0)], [c(0.0, 0.0), e]]
}
fn gate_sx() -> U2 {
    // SX = √X = (1/2)[[1+i, 1-i], [1-i, 1+i]]
    [[c(0.5,  0.5), c(0.5, -0.5)],
     [c(0.5, -0.5), c(0.5,  0.5)]]
}
fn gate_sxdg() -> U2 {
    // SX† = (√X)† = (1/2)[[1-i, 1+i], [1+i, 1-i]]
    [[c(0.5, -0.5), c(0.5,  0.5)],
     [c(0.5,  0.5), c(0.5, -0.5)]]
}

fn gate_rx(theta: f64) -> U2 {
    let cos = c((theta / 2.0).cos(), 0.0);
    let isin = c(0.0, -(theta / 2.0).sin());
    [[cos, isin], [isin, cos]]
}
fn gate_ry(theta: f64) -> U2 {
    let cos = c((theta / 2.0).cos(), 0.0);
    let sin = c((theta / 2.0).sin(), 0.0);
    [[cos, -sin], [sin, cos]]
}
fn gate_rz(theta: f64) -> U2 {
    let e_neg = Complex64::from_polar(1.0, -theta / 2.0);
    let e_pos = Complex64::from_polar(1.0, theta / 2.0);
    [[e_neg, c(0.0, 0.0)], [c(0.0, 0.0), e_pos]]
}
fn gate_phase(theta: f64) -> U2 {
    let e = Complex64::from_polar(1.0, theta);
    [[c(1.0, 0.0), c(0.0, 0.0)], [c(0.0, 0.0), e]]
}
fn gate_u(theta: f64, phi: f64, lambda: f64) -> U2 {
    let cos = (theta / 2.0).cos();
    let sin = (theta / 2.0).sin();
    [
        [c(cos, 0.0), -Complex64::from_polar(sin, lambda)],
        [Complex64::from_polar(sin, phi), Complex64::from_polar(cos, phi + lambda)],
    ]
}
fn gate_cu(theta: f64, phi: f64, lambda: f64, gamma: f64) -> U2 {
    let phase = Complex64::from_polar(1.0, gamma);
    let u = gate_u(theta, phi, lambda);
    [[phase * u[0][0], phase * u[0][1]], [phase * u[1][0], phase * u[1][1]]]
}

#[cfg(test)]
mod tests {
    use super::*;
    use cforge_core::{GateKind, Operation};

    fn bell_circuit() -> Circuit {
        let mut c = Circuit::new(2);
        c.push(Operation::new(GateKind::H, vec![0], vec![]));
        c.push(Operation::new(GateKind::Cx, vec![0, 1], vec![]));
        c
    }

    #[test]
    fn bell_state_amplitudes() {
        let backend = NativeStateVectorBackend;
        let result = backend.run(&bell_circuit(), 0, 0).unwrap();
        let sv = &result.statevector;
        let f = 1.0 / 2f64.sqrt();
        assert!((sv[0].re - f).abs() < 1e-10);
        assert!(sv[1].norm() < 1e-10);
        assert!(sv[2].norm() < 1e-10);
        assert!((sv[3].re - f).abs() < 1e-10);
    }

    #[test]
    fn x_gate_flips_zero_to_one() {
        let mut c = Circuit::new(1);
        c.push(Operation::new(GateKind::X, vec![0], vec![]));
        let result = NativeStateVectorBackend.run(&c, 0, 0).unwrap();
        let sv = &result.statevector;
        assert!(sv[0].norm() < 1e-10);
        assert!((sv[1].re - 1.0).abs() < 1e-10);
    }

    #[test]
    fn rz_returns_correct_phase() {
        let angle = std::f64::consts::FRAC_PI_2;
        let mut c = Circuit::new(1);
        c.push(Operation::new(GateKind::Rz, vec![0], vec![angle]));
        let result = NativeStateVectorBackend.run(&c, 0, 0).unwrap();
        let expected = Complex64::from_polar(1.0, -angle / 2.0);
        assert!((result.statevector[0] - expected).norm() < 1e-10);
    }

    #[test]
    fn shot_counts_bell_state() {
        let result = NativeStateVectorBackend.run(&bell_circuit(), 2000, 0).unwrap();
        let n00 = result.counts.get("00").copied().unwrap_or(0);
        let n11 = result.counts.get("11").copied().unwrap_or(0);
        assert_eq!(n00 + n11, 2000);
        // With 2000 shots the marginal error is ~2σ ≈ 45; allow 100.
        assert!((n00 as isize - 1000).unsigned_abs() < 200);
    }
}
