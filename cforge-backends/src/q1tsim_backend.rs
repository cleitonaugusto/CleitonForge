//! Simulation backend sourcing gate matrices from the `q1tsim` crate.
//!
//! Uses q1tsim's `Gate::matrix()` as the source of truth for unitary matrices
//! and applies them to a statevector using CleitonForge's own algorithm.
//!
//! **Rz convention**: q1tsim uses [[e^{-iλ/2},0],[0,e^{+iλ/2}]], which matches
//! CleitonForge's native statevector backend and roqoqo (IBM/Qiskit standard).
//! This is the third independent confirmation of the IBM convention — and the
//! third framework that disagrees with quantrs2-core on QAOA output states.
//!
//! **Gate coverage**: q1tsim ships controlled variants for almost every gate
//! (CH, CRX, CRY, CRZ, CU1, CU3, CCX, CV/CVdg) so fewer native fallbacks
//! are needed compared to the roqoqo backend.

use std::collections::HashMap;

use num_complex::Complex64;
use q1tsim::gates::{
    Gate, Sdg, Swap, Tdg, Vdg, CCX, CH, CRX, CRY, CRZ, CU1, CU3, CV, CX, CY, CZ, H, I, RX, RY, RZ,
    S, T, U1, U3, V, X, Y, Z,
};

use cforge_core::{Circuit, GateKind};

use crate::sample::sample_counts;
use crate::trait_def::{BackendError, SimulationBackend, SimulationResult};

const MAX_QUBITS: usize = 22;

/// Simulation backend that sources gate matrices from `q1tsim`.
pub struct Q1tSimBackend;

impl SimulationBackend for Q1tSimBackend {
    fn name(&self) -> &str {
        "statevector-q1tsim"
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
                "q1tsim backend supports up to {MAX_QUBITS} qubits, circuit has {n}"
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

// ── Matrix conversion helpers ─────────────────────────────────────────────────

type U2 = [[Complex64; 2]; 2];
type U4 = [[Complex64; 4]; 4];
type U8 = [[Complex64; 8]; 8];

// q1tsim uses num-complex 0.2; we use 0.4. The Complex types differ, so we
// extract .re/.im (both f64) and reconstruct our own Complex64.
#[inline]
fn q2c(z: &q1tsim::cmatrix::CNumber) -> Complex64 {
    Complex64::new(z.re, z.im)
}

fn m1(gate: impl Gate) -> Result<U2, String> {
    let m = gate.matrix();
    Ok([
        [q2c(&m[[0, 0]]), q2c(&m[[0, 1]])],
        [q2c(&m[[1, 0]]), q2c(&m[[1, 1]])],
    ])
}

fn m2(gate: impl Gate) -> Result<U4, String> {
    let m = gate.matrix();
    let mut out = [[Complex64::new(0.0, 0.0); 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            out[i][j] = q2c(&m[[i, j]]);
        }
    }
    Ok(out)
}

fn m3(gate: impl Gate) -> Result<U8, String> {
    let m = gate.matrix();
    let mut out = [[Complex64::new(0.0, 0.0); 8]; 8];
    for i in 0..8 {
        for j in 0..8 {
            out[i][j] = q2c(&m[[i, j]]);
        }
    }
    Ok(out)
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

/// q0 is more-significant bit in the 2-qubit subspace index (matches q1tsim convention).
fn apply2(sv: &mut [Complex64], q0: usize, q1: usize, u: U4) {
    let bit0 = 1 << q0;
    let bit1 = 1 << q1;
    for i in 0..sv.len() {
        if (i & bit0 != 0) || (i & bit1 != 0) {
            continue;
        }
        let idx = [i, i | bit1, i | bit0, i | bit0 | bit1];
        let old: [Complex64; 4] = [sv[idx[0]], sv[idx[1]], sv[idx[2]], sv[idx[3]]];
        for row in 0..4 {
            sv[idx[row]] = (0..4).map(|col| u[row][col] * old[col]).sum();
        }
    }
}

fn apply3(sv: &mut [Complex64], q0: usize, q1: usize, q2: usize, u: U8) {
    let bit0 = 1 << q0;
    let bit1 = 1 << q1;
    let bit2 = 1 << q2;
    for i in 0..sv.len() {
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

// ── Gate dispatch ─────────────────────────────────────────────────────────────

fn apply_gate(sv: &mut [Complex64], kind: GateKind, q: &[usize], p: &[f64]) -> Result<(), String> {
    match kind {
        GateKind::Id => apply1(sv, q[0], m1(I::new())?),
        GateKind::X => apply1(sv, q[0], m1(X::new())?),
        GateKind::Y => apply1(sv, q[0], m1(Y::new())?),
        GateKind::Z => apply1(sv, q[0], m1(Z::new())?),
        GateKind::H => apply1(sv, q[0], m1(H::new())?),
        GateKind::S => apply1(sv, q[0], m1(S::new())?),
        GateKind::Sdg => apply1(sv, q[0], m1(Sdg::new())?),
        GateKind::T => apply1(sv, q[0], m1(T::new())?),
        GateKind::Tdg => apply1(sv, q[0], m1(Tdg::new())?),
        GateKind::Sx => apply1(sv, q[0], m1(V::new())?),
        GateKind::Sxdg => apply1(sv, q[0], m1(Vdg::new())?),
        GateKind::Rx => apply1(sv, q[0], m1(RX::new(p[0]))?),
        GateKind::Ry => apply1(sv, q[0], m1(RY::new(p[0]))?),
        GateKind::Rz => apply1(sv, q[0], m1(RZ::new(p[0]))?),
        GateKind::Phase => apply1(sv, q[0], m1(U1::new(p[0]))?),
        GateKind::U => apply1(sv, q[0], m1(U3::new(p[0], p[1], p[2]))?),

        GateKind::Cx => apply2(sv, q[0], q[1], m2(CX::new())?),
        GateKind::Cy => apply2(sv, q[0], q[1], m2(CY::new())?),
        GateKind::Cz => apply2(sv, q[0], q[1], m2(CZ::new())?),
        GateKind::Ch => apply2(sv, q[0], q[1], m2(CH::new())?),
        GateKind::Csx => apply2(sv, q[0], q[1], m2(CV::new())?),
        GateKind::Crx => apply2(sv, q[0], q[1], m2(CRX::new(p[0]))?),
        GateKind::Cry => apply2(sv, q[0], q[1], m2(CRY::new(p[0]))?),
        GateKind::Crz => apply2(sv, q[0], q[1], m2(CRZ::new(p[0]))?),
        GateKind::Cp => apply2(sv, q[0], q[1], m2(CU1::new(p[0]))?),
        GateKind::Swap => apply2(sv, q[0], q[1], m2(Swap::new())?),

        // CU: q1tsim's CU3 matches our Cu gate (theta, phi, lambda + global phase γ ignored)
        GateKind::Cu => {
            let u = m2(CU3::new(p[0], p[1], p[2]))?;
            apply2(sv, q[0], q[1], u);
        }

        GateKind::Ccx => apply3(sv, q[0], q[1], q[2], m3(CCX::new())?),
        GateKind::Cswap => apply_cswap(sv, q[0], q[1], q[2]),
    }
    Ok(())
}

fn apply_cswap(sv: &mut [Complex64], ctrl: usize, a: usize, b: usize) {
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
    use crate::{NativeStateVectorBackend, SimulationBackend};
    use cforge_core::Operation;

    fn fidelity(a: &[Complex64], b: &[Complex64]) -> f64 {
        let inner: Complex64 = a.iter().zip(b).map(|(x, y)| x.conj() * y).sum();
        inner.norm_sqr()
    }

    fn run_both(circuit: &Circuit) -> (Vec<Complex64>, Vec<Complex64>) {
        let r1 = NativeStateVectorBackend.run(circuit, 0, 0).unwrap();
        let r2 = Q1tSimBackend.run(circuit, 0, 0).unwrap();
        (r1.statevector, r2.statevector)
    }

    #[test]
    fn bell_state_q1tsim_agrees() {
        let mut c = Circuit::new(2);
        c.push(Operation::new(GateKind::H, vec![0], vec![]));
        c.push(Operation::new(GateKind::Cx, vec![0, 1], vec![]));
        let (sv1, sv2) = run_both(&c);
        assert!(fidelity(&sv1, &sv2) > 0.9999);
    }

    #[test]
    fn rz_on_superposition_q1tsim_agrees_with_native() {
        // q1tsim Rz = [[e^{-iθ/2},0],[0,e^{+iθ/2}]] — same IBM convention.
        // Third confirmation: native ≡ roqoqo ≡ q1tsim ≠ quantrs2.
        let mut c = Circuit::new(1);
        c.push(Operation::new(GateKind::H, vec![0], vec![]));
        c.push(Operation::new(
            GateKind::Rz,
            vec![0],
            vec![std::f64::consts::FRAC_PI_4],
        ));
        let (sv1, sv2) = run_both(&c);
        assert!(
            fidelity(&sv1, &sv2) > 0.9999,
            "fidelity = {}",
            fidelity(&sv1, &sv2)
        );
    }

    #[test]
    fn qaoa_maxcut_q1tsim_agrees_with_native() {
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
        assert!(f > 0.9999, "QAOA q1tsim vs native fidelity = {f}");
    }

    #[test]
    fn toffoli_q1tsim_agrees() {
        let mut c = Circuit::new(3);
        for q in 0..3 {
            c.push(Operation::new(GateKind::H, vec![q], vec![]));
        }
        c.push(Operation::new(GateKind::Ccx, vec![0, 1, 2], vec![]));
        let (sv1, sv2) = run_both(&c);
        assert!(fidelity(&sv1, &sv2) > 0.9999);
    }

    #[test]
    fn q1tsim_backend_name() {
        assert_eq!(Q1tSimBackend.name(), "statevector-q1tsim");
    }

    #[test]
    fn q1tsim_all_controlled_gates() {
        // Smoke test: all controlled gates don't panic or error
        let mut c = Circuit::new(3);
        c.push(Operation::new(GateKind::Ch, vec![0, 1], vec![]));
        c.push(Operation::new(GateKind::Csx, vec![0, 1], vec![]));
        c.push(Operation::new(GateKind::Crx, vec![0, 1], vec![0.5]));
        c.push(Operation::new(GateKind::Cry, vec![0, 1], vec![0.5]));
        c.push(Operation::new(GateKind::Crz, vec![0, 1], vec![0.5]));
        c.push(Operation::new(GateKind::Cp, vec![0, 1], vec![0.5]));
        c.push(Operation::new(GateKind::Swap, vec![1, 2], vec![]));
        let r = Q1tSimBackend.run(&c, 0, 0);
        assert!(r.is_ok());
    }
}
