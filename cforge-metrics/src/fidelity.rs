//! Fidelity computation between simulation results and reference states.

use num_complex::Complex64;

/// Computes the fidelity between two pure statevectors: |⟨ψ_ref|ψ⟩|².
///
/// Returns `None` when the vectors have different lengths (incompatible
/// circuits) or when either has zero norm (degenerate state).
pub fn statevector_fidelity(result: &[Complex64], reference: &[Complex64]) -> Option<f64> {
    if result.len() != reference.len() || result.is_empty() {
        return None;
    }
    let norm_result_sq: f64 = result.iter().map(|a| a.norm_sqr()).sum();
    let norm_ref_sq: f64 = reference.iter().map(|a| a.norm_sqr()).sum();
    if norm_result_sq < 1e-15 || norm_ref_sq < 1e-15 {
        return None;
    }
    let inner: Complex64 = reference
        .iter()
        .zip(result)
        .map(|(r, s)| r.conj() * s)
        .sum();
    Some(inner.norm_sqr() / (norm_result_sq * norm_ref_sq))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cx(re: f64, im: f64) -> Complex64 {
        Complex64::new(re, im)
    }

    #[test]
    fn identical_states_have_fidelity_one() {
        let f = 1.0 / 2f64.sqrt();
        let sv = vec![cx(f, 0.0), cx(0.0, 0.0), cx(0.0, 0.0), cx(f, 0.0)];
        assert!((statevector_fidelity(&sv, &sv).unwrap() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn orthogonal_states_have_fidelity_zero() {
        let zero = vec![cx(1.0, 0.0), cx(0.0, 0.0)];
        let one = vec![cx(0.0, 0.0), cx(1.0, 0.0)];
        assert!(statevector_fidelity(&zero, &one).unwrap() < 1e-15);
    }

    #[test]
    fn different_lengths_return_none() {
        let a = vec![cx(1.0, 0.0)];
        let b = vec![cx(1.0, 0.0), cx(0.0, 0.0)];
        assert!(statevector_fidelity(&a, &b).is_none());
    }
}
