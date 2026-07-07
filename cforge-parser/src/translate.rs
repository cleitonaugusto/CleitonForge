//! Mapping from gate name strings to cforge-core `GateKind`,
//! including parameter rewriting for aliased gates.

use cforge_core::GateKind;
use std::f64::consts::FRAC_PI_2;

/// Translates a gate name and its already-evaluated parameters into the
/// canonical `(GateKind, params)` pair used in the circuit IR.
///
/// Handles:
/// - Direct name→kind mappings (case-insensitive, QASM2 + QASM3 aliases)
/// - Parameter rewriting for compound aliases:
///   - `u2(φ, λ)` → `U(π/2, φ, λ)`
///
/// Returns `None` for gate names that are not part of the supported gate set.
pub fn translate_gate(name: &str, params: Vec<f64>) -> Option<(GateKind, Vec<f64>)> {
    match name.to_lowercase().as_str() {
        // u2(phi, lam) is shorthand for u(pi/2, phi, lam).
        // IBM/Qiskit circuits use this heavily; rewrite at parse time so
        // the IR stays clean (no U2 variant needed in GateKind).
        "u2" => {
            let mut full = Vec::with_capacity(3);
            full.push(FRAC_PI_2);
            full.extend(params);
            Some((GateKind::U, full))
        }
        other => gate_kind_from_name(other).map(|kind| (kind, params)),
    }
}

/// Looks up a gate name (case-insensitive) and returns the matching `GateKind`.
///
/// Does not handle parameter rewriting; use [`translate_gate`] when you
/// also have evaluated parameters available.
pub fn gate_kind_from_name(name: &str) -> Option<GateKind> {
    match name.to_lowercase().as_str() {
        "id" | "i" => Some(GateKind::Id),
        "x" => Some(GateKind::X),
        "y" => Some(GateKind::Y),
        "z" => Some(GateKind::Z),
        "h" => Some(GateKind::H),
        "s" => Some(GateKind::S),
        "sdg" => Some(GateKind::Sdg),
        "t" => Some(GateKind::T),
        "tdg" => Some(GateKind::Tdg),
        "sx" => Some(GateKind::Sx),
        "sxdg" => Some(GateKind::Sxdg),
        "rx" => Some(GateKind::Rx),
        "ry" => Some(GateKind::Ry),
        "rz" => Some(GateKind::Rz),
        "p" | "phase" | "u1" => Some(GateKind::Phase),
        "u" | "u3" => Some(GateKind::U),
        "cx" | "cnot" => Some(GateKind::Cx),
        "cy" => Some(GateKind::Cy),
        "cz" => Some(GateKind::Cz),
        "ch" => Some(GateKind::Ch),
        "csx" => Some(GateKind::Csx),
        "crx" => Some(GateKind::Crx),
        "cry" => Some(GateKind::Cry),
        "crz" => Some(GateKind::Crz),
        "cp" | "cphase" => Some(GateKind::Cp),
        "cu" | "cu3" => Some(GateKind::Cu),
        "swap" => Some(GateKind::Swap),
        "ccx" | "toffoli" | "ccnot" => Some(GateKind::Ccx),
        "cswap" | "fredkin" => Some(GateKind::Cswap),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn u2_rewrites_to_u_with_pi_over_2() {
        let (kind, params) = translate_gate("u2", vec![0.5, 1.0]).unwrap();
        assert_eq!(kind, GateKind::U);
        assert_eq!(params.len(), 3);
        assert!((params[0] - FRAC_PI_2).abs() < 1e-15);
        assert_eq!(params[1], 0.5);
        assert_eq!(params[2], 1.0);
    }

    #[test]
    fn u2_uppercase_also_works() {
        let (kind, _) = translate_gate("U2", vec![0.0, 0.0]).unwrap();
        assert_eq!(kind, GateKind::U);
    }

    #[test]
    fn direct_gates_pass_params_unchanged() {
        let (kind, params) = translate_gate("rx", vec![1.5]).unwrap();
        assert_eq!(kind, GateKind::Rx);
        assert_eq!(params, vec![1.5]);
    }

    #[test]
    fn unknown_gate_returns_none() {
        assert!(translate_gate("fakegate", vec![]).is_none());
    }
}
