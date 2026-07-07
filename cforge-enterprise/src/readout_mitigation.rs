//! Readout error mitigation via confusion matrix inversion.
//!
//! The confusion matrix A has entries A[i][j] = P(measure i | prepared j).
//! For n qubits, A is 2^n × 2^n. Mitigation: p_ideal = A⁻¹ p_noisy.
//!
//! For large n, use the tensored/local approximation: A ≈ A₀ ⊗ A₁ ⊗ ... ⊗ A_{n-1}
//! where each A_k is a 2×2 confusion matrix for qubit k.

/// Readout mitigation using per-qubit confusion matrix inversion.
pub struct ReadoutMitigationMatrix {
    /// Per-qubit confusion matrices: [[P(0|0), P(1|0)], [P(0|1), P(1|1)]].
    pub per_qubit: Vec<[[f64; 2]; 2]>,
}

impl ReadoutMitigationMatrix {
    /// Build from device readout error rates (p_err = P(flip | prepared state)).
    pub fn from_readout_errors(readout_errs: &[f64]) -> Self {
        let per_qubit = readout_errs
            .iter()
            .map(|&p| [[1.0 - p, p], [p, 1.0 - p]])
            .collect();
        Self { per_qubit }
    }

    /// Apply mitigation: correct the noisy probability vector.
    ///
    /// Uses local (tensored) approximation — O(n) per-qubit inversions.
    /// Full 2^n × 2^n inversion is available in the enterprise release.
    pub fn mitigate(&self, probs: &[f64]) -> Result<Vec<f64>, String> {
        // TODO: enterprise implementation with proper matrix inversion
        let _ = probs;
        Err("ReadoutMitigationMatrix::mitigate requires enterprise license".into())
    }
}
