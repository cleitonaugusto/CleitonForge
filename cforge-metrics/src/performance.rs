//! Wall-clock timing and peak-memory wrappers around backend calls.

use std::time::Instant;
use cforge_backends::{BackendError, SimulationBackend, SimulationResult};
use cforge_core::{Circuit, MetricsResult};

use crate::circuit_stats::compute_stats;
use crate::fidelity::statevector_fidelity;

/// Runs `backend` on `circuit`, measures wall-clock time, and produces
/// a fully-populated [`MetricsResult`].
///
/// `reference` is the expected statevector; pass `None` when the
/// expected output is not known (fidelity will be `None`).
pub fn measure(
    backend: &dyn SimulationBackend,
    circuit: &Circuit,
    shots: usize,
    reference: Option<&[num_complex::Complex64]>,
) -> Result<MetricsResult, BackendError> {
    let stats = compute_stats(circuit);

    let t0 = Instant::now();
    let result: SimulationResult = backend.run(circuit, shots)?;
    let elapsed_ms = t0.elapsed().as_secs_f64() * 1000.0;

    let fidelity = reference.and_then(|r| statevector_fidelity(&result.statevector, r));

    Ok(MetricsResult {
        fidelity,
        depth: stats.depth,
        gate_count: stats.gate_count,
        execution_time_ms: elapsed_ms,
        memory_bytes: None, // platform-independent peak-memory requires OS-specific APIs
        backend_name: backend.name().to_string(),
    })
}
