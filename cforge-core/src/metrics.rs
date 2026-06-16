//! Standardized result types for circuit and backend measurements. The
//! actual computation logic lives in `cforge-metrics`; this crate only
//! defines the shared shape of the result so every backend reports
//! comparable data.

/// Standardized metrics for a single circuit run on a single backend.
#[derive(Debug, Clone, PartialEq)]
pub struct MetricsResult {
    /// Fidelity against a reference state, when one is known.
    pub fidelity: Option<f64>,
    /// Circuit depth (number of sequential layers of gates), computed on
    /// the canonical IR so it's comparable across backends.
    pub depth: usize,
    /// Total number of gates in the circuit, computed on the canonical IR.
    pub gate_count: usize,
    pub execution_time_ms: f64,
    /// Peak memory used during the backend call, when measurable.
    pub memory_bytes: Option<u64>,
    pub backend_name: String,
}
