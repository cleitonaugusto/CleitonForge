//! The `SimulationBackend` trait and associated result/error types.

use std::collections::HashMap;
use num_complex::Complex64;
use cforge_core::Circuit;

/// Measurement outcome and statevector for one circuit run.
#[derive(Debug, Clone)]
pub struct SimulationResult {
    /// Final amplitude vector (|ψ⟩) after all gates have been applied.
    pub statevector: Vec<Complex64>,
    /// Bitstring measurement counts, populated when `shots > 0`.
    pub counts: HashMap<String, usize>,
}

impl SimulationResult {
    /// Returns the probability for each computational basis state.
    pub fn probabilities(&self) -> Vec<f64> {
        self.statevector.iter().map(|a| a.norm_sqr()).collect()
    }

    /// Computes the fidelity between this result's statevector and a
    /// reference statevector `|ref⟩`: |⟨ref|ψ⟩|².
    pub fn fidelity(&self, reference: &[Complex64]) -> Option<f64> {
        if self.statevector.len() != reference.len() {
            return None;
        }
        let inner: Complex64 = reference
            .iter()
            .zip(self.statevector.iter())
            .map(|(r, s)| r.conj() * s)
            .sum();
        Some(inner.norm_sqr())
    }
}

/// An error from a simulation backend.
#[derive(Debug, Clone)]
pub struct BackendError(pub String);

impl std::fmt::Display for BackendError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "backend error: {}", self.0)
    }
}

impl std::error::Error for BackendError {}

/// A pluggable quantum simulation backend.
///
/// Implementations must be stateless — all state lives in the `Circuit`
/// and the returned `SimulationResult`. This allows multiple backends to
/// be benchmarked on the identical input without interference.
pub trait SimulationBackend {
    /// A short, stable name identifying this backend (e.g. `"statevector-native"`).
    fn name(&self) -> &str;

    /// Simulate `circuit` and return the final quantum state.
    ///
    /// When `shots > 0`, also sample `shots` measurements from the
    /// resulting probability distribution and populate `counts`.
    fn run(
        &self,
        circuit: &Circuit,
        shots: usize,
    ) -> Result<SimulationResult, BackendError>;
}
