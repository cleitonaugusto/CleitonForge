//! Automatic classification of a minimized divergence.
//!
//! Triage answers the two questions a maintainer asks first:
//! *what kind of bug is this* (which QGCS dimension) and *could any
//! standard benchmark have caught it* (benchmark visibility). The
//! second answer comes directly from the oracle hierarchy: a divergence
//! at N1-amplitude with agreement at N2-probability is invisible to
//! every sampling benchmark by the conjugation-invariance theorem.

use cforge_core::{Circuit, GateKind};

use crate::oracle::{amplitude_distance_mod_phase, probability_distance};
use num_complex::Complex64;

/// Coarse bug classification aligned with QGCS dimensions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BugClass {
    /// Divergence in the diagonal rotation family (Rz/Phase/Crz/Cp) —
    /// rotation sign or phase-placement convention.
    RotationConvention,
    /// Divergence involving discrete phase gates (S/T/Sx and daggers).
    DiscretePhaseGate,
    /// Divergence only visible with multi-qubit gates — entangler or
    /// qubit-ordering (endianness) convention.
    MultiQubitConvention,
    /// Anything else.
    Unclassified,
}

impl BugClass {
    pub fn label(&self) -> &'static str {
        match self {
            BugClass::RotationConvention => "rotation-convention (QGCS dim: Rz/Phase sign)",
            BugClass::DiscretePhaseGate => "discrete-phase-gate (QGCS dim: S/T/√X phase)",
            BugClass::MultiQubitConvention => "multi-qubit convention (QGCS dim: ordering)",
            BugClass::Unclassified => "unclassified",
        }
    }
}

/// Full triage verdict for a minimized counterexample.
#[derive(Debug, Clone)]
pub struct Triage {
    pub class: BugClass,
    /// N1 distance (amplitude, mod global phase).
    pub amplitude_distance: f64,
    /// N2 distance (probability).
    pub probability_distance: f64,
    /// True when some sampling benchmark could in principle observe the
    /// divergence (N2 > tol). False = invisible to QV/XEB/mirror/RB/QPE.
    pub sampling_visible: bool,
}

impl Triage {
    pub fn visibility_label(&self) -> &'static str {
        if self.sampling_visible {
            "visible to sampling benchmarks (probability-level)"
        } else {
            "INVISIBLE to all sampling benchmarks (amplitude-only — theorem class)"
        }
    }
}

/// Classifies a minimal witness given the two backends' statevectors.
pub fn triage(
    minimal: &Circuit,
    sv_reference: &[Complex64],
    sv_device: &[Complex64],
    tol: f64,
) -> Triage {
    let amp = amplitude_distance_mod_phase(sv_reference, sv_device);
    let prob = probability_distance(sv_reference, sv_device);

    let has_rotation = minimal.operations.iter().any(|op| {
        matches!(
            op.kind,
            GateKind::Rz | GateKind::Phase | GateKind::Crz | GateKind::Cp
        )
    });
    let has_discrete_phase = minimal.operations.iter().any(|op| {
        matches!(
            op.kind,
            GateKind::S
                | GateKind::Sdg
                | GateKind::T
                | GateKind::Tdg
                | GateKind::Sx
                | GateKind::Sxdg
                | GateKind::Csx
        )
    });
    let only_multi_qubit_complex = minimal
        .operations
        .iter()
        .all(|op| op.kind.num_qubits() >= 2 || op.kind == GateKind::H);

    let class = if has_rotation {
        BugClass::RotationConvention
    } else if has_discrete_phase {
        BugClass::DiscretePhaseGate
    } else if only_multi_qubit_complex {
        BugClass::MultiQubitConvention
    } else {
        BugClass::Unclassified
    };

    Triage {
        class,
        amplitude_distance: amp,
        probability_distance: prob,
        sampling_visible: prob > tol,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cforge_core::Operation;

    #[test]
    fn conjugated_state_is_flagged_benchmark_invisible() {
        let f = std::f64::consts::FRAC_1_SQRT_2;
        let a = vec![Complex64::new(f, 0.0), Complex64::new(0.0, f)];
        let b = vec![Complex64::new(f, 0.0), Complex64::new(0.0, -f)];
        let mut c = Circuit::new(1);
        c.push(Operation::new(GateKind::H, vec![0], vec![]));
        c.push(Operation::new(GateKind::Rz, vec![0], vec![0.5]));

        let t = triage(&c, &a, &b, 1e-6);
        assert!(!t.sampling_visible);
        assert_eq!(t.class, BugClass::RotationConvention);
    }
}
