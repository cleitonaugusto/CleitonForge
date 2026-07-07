use cforge_core::{Circuit, GateKind, Operation};

/// The sign convention used for Rz-family gates.
///
/// The IBM/Qiskit standard defines `Rz(λ) = [[e^{-iλ/2}, 0], [0, e^{+iλ/2}]]`.
/// Some frameworks (notably quantrs2-core) use the opposite sign.
///
/// Convention divergence is invisible in Clifford-only or random-angle circuits
/// but produces zero cross-backend fidelity on circuits with specific angles
/// (QAOA, VQE). Empirically confirmed by CleitonForge benchmark suite
/// across 4 statevector backends.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RzConvention {
    /// IBM / Qiskit / OpenQASM 3 standard.
    ///
    /// `Rz(λ) = [[e^{-iλ/2}, 0], [0, e^{+iλ/2}]]`
    ///
    /// Empirically confirmed: CleitonForge native, roqoqo, q1tsim.
    Standard,

    /// Opposite sign convention.
    ///
    /// `Rz(λ) = [[e^{+iλ/2}, 0], [0, e^{-iλ/2}]]`
    ///
    /// Empirically confirmed: quantrs2-core.
    Reversed,
}

impl RzConvention {
    /// Returns the convention name as a static string.
    pub fn name(self) -> &'static str {
        match self {
            Self::Standard => "standard (IBM/Qiskit)",
            Self::Reversed => "reversed (quantrs2)",
        }
    }
}

/// Rewrites `circuit` so that its Rz-family gate angles match `to`.
///
/// When `from == to` the original circuit is returned unchanged (cloned).
/// The transformation is its own inverse: applying it twice yields the
/// original circuit.
///
/// # Gate coverage
///
/// | Gate | Params negated | Basis |
/// |------|----------------|-------|
/// | `Rz(λ)` | `λ` | empirically confirmed |
/// | `Phase(λ)` | `λ` | `P = e^{iλ/2} Rz(λ)` up to global phase |
/// | `Crz(λ)` | `λ` | controlled `Rz` |
/// | `Cp(λ)` | `λ` | controlled `Phase` |
/// | `U(θ,φ,λ)` | `φ`, `λ` | `θ` is `Ry`-style, unaffected |
/// | `Cu(θ,φ,λ,γ)` | `φ`, `λ`, `γ` | controlled `U` with global phase |
///
/// `Rx`, `Ry`, `H`, `Cx`, and all Clifford gates are left unchanged.
pub fn normalize_convention(circuit: &Circuit, from: RzConvention, to: RzConvention) -> Circuit {
    if from == to {
        return circuit.clone();
    }
    let mut out = Circuit::new(circuit.num_qubits());
    for op in &circuit.operations {
        out.push(normalized_op(op));
    }
    out
}

fn normalized_op(op: &Operation) -> Operation {
    use GateKind::*;
    let mut out = op.clone();
    match op.kind {
        Rz | Phase | Crz | Cp => {
            out.params[0] = -op.params[0];
        }
        U => {
            // params: [θ, φ, λ] — θ is Ry-style (unaffected)
            out.params[1] = -op.params[1]; // φ
            out.params[2] = -op.params[2]; // λ
        }
        Cu => {
            // params: [θ, φ, λ, γ] — θ is Ry-style (unaffected)
            out.params[1] = -op.params[1]; // φ
            out.params[2] = -op.params[2]; // λ
            out.params[3] = -op.params[3]; // γ (global phase offset)
        }
        _ => {}
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use cforge_core::Operation;
    use std::f64::consts::PI;

    fn rz_circuit(lambda: f64) -> Circuit {
        let mut c = Circuit::new(1);
        c.push(Operation::new(GateKind::Rz, vec![0], vec![lambda]));
        c
    }

    #[test]
    fn noop_when_same_convention() {
        let original = rz_circuit(PI / 4.0);
        let normalized =
            normalize_convention(&original, RzConvention::Standard, RzConvention::Standard);
        assert_eq!(
            original.operations[0].params[0],
            normalized.operations[0].params[0]
        );
    }

    #[test]
    fn rz_angle_negated_on_convention_flip() {
        let circuit = rz_circuit(PI / 4.0);
        let normalized =
            normalize_convention(&circuit, RzConvention::Reversed, RzConvention::Standard);
        assert!((normalized.operations[0].params[0] - (-PI / 4.0)).abs() < 1e-15);
    }

    #[test]
    fn double_normalization_is_identity() {
        let circuit = rz_circuit(PI / 3.0);
        let once = normalize_convention(&circuit, RzConvention::Reversed, RzConvention::Standard);
        let twice = normalize_convention(&once, RzConvention::Reversed, RzConvention::Standard);
        assert!((twice.operations[0].params[0] - circuit.operations[0].params[0]).abs() < 1e-15);
    }

    #[test]
    fn phase_gate_angle_negated() {
        let mut c = Circuit::new(1);
        c.push(Operation::new(GateKind::Phase, vec![0], vec![PI / 2.0]));
        let n = normalize_convention(&c, RzConvention::Reversed, RzConvention::Standard);
        assert!((n.operations[0].params[0] - (-PI / 2.0)).abs() < 1e-15);
    }

    #[test]
    fn u_gate_phi_lambda_negated_theta_preserved() {
        let theta = 1.1;
        let phi = 0.7;
        let lambda = -0.3;
        let mut c = Circuit::new(1);
        c.push(Operation::new(
            GateKind::U,
            vec![0],
            vec![theta, phi, lambda],
        ));
        let n = normalize_convention(&c, RzConvention::Reversed, RzConvention::Standard);
        let p = &n.operations[0].params;
        assert!((p[0] - theta).abs() < 1e-15); // θ unchanged
        assert!((p[1] - (-phi)).abs() < 1e-15);
        assert!((p[2] - (-lambda)).abs() < 1e-15);
    }

    #[test]
    fn non_rz_gates_unchanged() {
        let mut c = Circuit::new(2);
        c.push(Operation::new(GateKind::H, vec![0], vec![]));
        c.push(Operation::new(GateKind::Cx, vec![0, 1], vec![]));
        c.push(Operation::new(GateKind::Rx, vec![0], vec![PI / 2.0]));
        c.push(Operation::new(GateKind::Ry, vec![1], vec![PI / 3.0]));
        let n = normalize_convention(&c, RzConvention::Reversed, RzConvention::Standard);
        assert_eq!(c.operations, n.operations);
    }

    #[test]
    fn crz_and_cp_negated() {
        let mut c = Circuit::new(2);
        c.push(Operation::new(GateKind::Crz, vec![0, 1], vec![PI]));
        c.push(Operation::new(GateKind::Cp, vec![0, 1], vec![PI / 4.0]));
        let n = normalize_convention(&c, RzConvention::Reversed, RzConvention::Standard);
        assert!((n.operations[0].params[0] - (-PI)).abs() < 1e-15);
        assert!((n.operations[1].params[0] - (-PI / 4.0)).abs() < 1e-15);
    }
}
