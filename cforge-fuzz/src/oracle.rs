//! The three-level differential oracle hierarchy.
//!
//! Each level compares two final statevectors with strictly less
//! information than the level below it:
//!
//! - **N1 Amplitude** — equality of amplitudes modulo one global phase.
//!   The strongest observable-in-principle notion of "same circuit
//!   semantics".
//! - **N2 Probability** — equality of |amplitude|² only. This is
//!   everything any sampling benchmark (QV, XEB, mirror, RB, QPE
//!   histograms) can ever see.
//! - **N3 Observable** — equality of single-qubit ⟨Z⟩ expectations, the
//!   coarsest practical signal.
//!
//! The hierarchy is strict: a consistent conjugation bug diverges at N1
//! while being *provably* identical at N2 and N3 (conjugation-invariance
//! theorem). A divergence at N1 with agreement at N2 is therefore the
//! fingerprint of a bug class invisible to the entire standard
//! benchmark stack.

use num_complex::Complex64;

/// Oracle strictness level.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OracleLevel {
    /// N1 — amplitude equality modulo global phase.
    Amplitude,
    /// N2 — probability (|amplitude|²) equality.
    Probability,
    /// N3 — per-qubit ⟨Z⟩ expectation equality.
    Observable,
}

impl OracleLevel {
    pub fn label(&self) -> &'static str {
        match self {
            OracleLevel::Amplitude => "N1-amplitude",
            OracleLevel::Probability => "N2-probability",
            OracleLevel::Observable => "N3-observable",
        }
    }
}

/// Distance between two statevectors at the given oracle level.
/// Returns `f64::INFINITY` on dimension mismatch.
pub fn distance(a: &[Complex64], b: &[Complex64], level: OracleLevel) -> f64 {
    if a.len() != b.len() || a.is_empty() {
        return f64::INFINITY;
    }
    match level {
        OracleLevel::Amplitude => amplitude_distance_mod_phase(a, b),
        OracleLevel::Probability => probability_distance(a, b),
        OracleLevel::Observable => observable_z_distance(a, b),
    }
}

/// N1: max_i |e^{iφ}·a_i − b_i| minimized over the global phase φ.
///
/// The optimal alignment is φ = arg⟨a|b⟩, which maximizes
/// Re(e^{iφ}⟨a|b⟩); a legitimate implementation difference of a pure
/// global phase then yields distance 0, so it is never reported as a
/// bug.
pub fn amplitude_distance_mod_phase(a: &[Complex64], b: &[Complex64]) -> f64 {
    let overlap: Complex64 = a.iter().zip(b).map(|(x, y)| x.conj() * y).sum();
    let phase = if overlap.norm() > 1e-15 {
        overlap / overlap.norm()
    } else {
        // Orthogonal states: no phase can align them; any choice
        // exposes the (maximal) divergence.
        Complex64::new(1.0, 0.0)
    };
    a.iter()
        .zip(b)
        .map(|(x, y)| (x * phase - y).norm())
        .fold(0.0_f64, f64::max)
}

/// N2: max_i ||a_i|² − |b_i|²| — everything sampling can see.
pub fn probability_distance(a: &[Complex64], b: &[Complex64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| (x.norm_sqr() - y.norm_sqr()).abs())
        .fold(0.0_f64, f64::max)
}

/// N3: max over qubits of |⟨Z_q⟩_a − ⟨Z_q⟩_b|.
pub fn observable_z_distance(a: &[Complex64], b: &[Complex64]) -> f64 {
    let n = a.len().trailing_zeros() as usize;
    (0..n)
        .map(|q| (expect_z(a, q) - expect_z(b, q)).abs())
        .fold(0.0_f64, f64::max)
}

fn expect_z(sv: &[Complex64], qubit: usize) -> f64 {
    let bit = 1usize << qubit;
    sv.iter()
        .enumerate()
        .map(|(i, amp)| {
            let sign = if i & bit == 0 { 1.0 } else { -1.0 };
            sign * amp.norm_sqr()
        })
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn c(re: f64, im: f64) -> Complex64 {
        Complex64::new(re, im)
    }

    #[test]
    fn hierarchy_is_strict_on_conjugated_state() {
        // |ψ⟩ = (|0⟩ + i|1⟩)/√2 vs its conjugate (|0⟩ − i|1⟩)/√2:
        // different at N1, identical at N2 and N3.
        let f = std::f64::consts::FRAC_1_SQRT_2;
        let a = vec![c(f, 0.0), c(0.0, f)];
        let b = vec![c(f, 0.0), c(0.0, -f)];
        assert!(distance(&a, &b, OracleLevel::Amplitude) > 0.5);
        assert!(distance(&a, &b, OracleLevel::Probability) < 1e-15);
        assert!(distance(&a, &b, OracleLevel::Observable) < 1e-15);
    }

    #[test]
    fn identical_states_have_zero_distance_everywhere() {
        let a = vec![c(0.6, 0.0), c(0.0, 0.8)];
        for level in [
            OracleLevel::Amplitude,
            OracleLevel::Probability,
            OracleLevel::Observable,
        ] {
            assert!(distance(&a, &a, level) < 1e-15);
        }
    }

    #[test]
    fn dimension_mismatch_is_infinite() {
        let a = vec![c(1.0, 0.0)];
        let b = vec![c(1.0, 0.0), c(0.0, 0.0)];
        assert_eq!(distance(&a, &b, OracleLevel::Amplitude), f64::INFINITY);
    }
}
