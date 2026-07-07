//! Zero-Noise Extrapolation (ZNE) — error mitigation skeleton.
//!
//! ZNE artificially scales the noise level by a factor λ > 1 (by folding gates:
//! G → G G† G), measures at several noise scales, then extrapolates back to λ=0.
//!
//! This is one of the most practical NISQ error mitigation techniques — it
//! does not require ancilla qubits or post-selection.
//!
//! ## Reference
//!
//! Temme, K., Bravyi, S. & Gambetta, J.M. "Error mitigation for short-depth
//! quantum circuits." PRL 119, 180509 (2017).

use cforge_core::Circuit;

/// ZNE estimator: fold gates at several noise scales and extrapolate to zero.
pub struct ZneEstimator {
    /// Noise scale factors λ ≥ 1. Example: [1.0, 1.5, 2.0, 2.5, 3.0].
    pub scale_factors: Vec<f64>,
    /// Extrapolation method.
    pub method: ExtrapolationMethod,
}

#[derive(Debug, Clone, Copy)]
pub enum ExtrapolationMethod {
    Linear,
    Exponential,
    Richardson,
}

impl ZneEstimator {
    pub fn new(scale_factors: Vec<f64>) -> Self {
        Self {
            scale_factors,
            method: ExtrapolationMethod::Richardson,
        }
    }

    /// Fold a circuit to scale noise by factor λ (integer scale only in v1).
    ///
    /// Gate folding: replace every gate G with G (G† G)^{(λ-1)/2}.
    /// For integer odd λ: 1x, 3x, 5x folding.
    pub fn fold_circuit(&self, circuit: &Circuit, scale: f64) -> Circuit {
        // TODO: enterprise implementation
        // This skeleton documents the API surface for the enterprise release.
        let _ = scale;
        circuit.clone()
    }

    /// Run ZNE: simulate at each scale factor, extrapolate to zero noise.
    ///
    /// Returns the extrapolated expectation value of observable `observable_idx`
    /// (index into the probability distribution, treated as a Pauli Z eigenvalue
    /// ±1 measurement).
    pub fn estimate(
        &self,
        _circuit: &Circuit,
        _observable_idx: usize,
        _shots: usize,
    ) -> Result<f64, String> {
        Err("ZneEstimator::estimate requires enterprise license".into())
    }
}
