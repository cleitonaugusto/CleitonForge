//! Wall-clock timing and peak-memory measurement around backend calls.

use std::time::Instant;
use cforge_backends::{BackendError, SimulationBackend, SimulationResult};
use cforge_core::{Circuit, MetricsResult};

use crate::circuit_stats::compute_stats;
use crate::fidelity::statevector_fidelity;
use crate::memory;

/// Runs `backend` on `circuit`, measures wall-clock time and peak memory,
/// and produces a fully-populated [`MetricsResult`].
///
/// **Memory**: On Linux the RSS delta is read from `/proc/self/status`
/// while the statevector is still live in memory, giving the actual OS
/// footprint. On other platforms (or for circuits smaller than one OS
/// page) the theoretical minimum `2^n × 16 bytes` is used instead.
///
/// `reference` is the expected statevector for fidelity; pass `None`
/// when not known.
pub fn measure(
    backend: &dyn SimulationBackend,
    circuit: &Circuit,
    shots: usize,
    seed: u64,
    reference: Option<&[num_complex::Complex64]>,
) -> Result<MetricsResult, BackendError> {
    let stats = compute_stats(circuit);
    let theoretical = memory::statevector_memory_bytes(circuit.num_qubits());

    // Sample RSS before the allocation.
    let rss_before = memory::current_rss_bytes();

    let t0 = Instant::now();
    let result: SimulationResult = backend.run(circuit, shots, seed)?;
    let elapsed_ms = t0.elapsed().as_secs_f64() * 1000.0;

    // Sample RSS while the statevector is still live inside `result`.
    let rss_after = memory::current_rss_bytes();

    let fidelity = reference.and_then(|r| statevector_fidelity(&result.statevector, r));

    // Use RSS delta when the allocation is large enough to register in
    // page-granular OS accounting; otherwise the theoretical minimum is
    // more accurate than a noisy sub-page delta.
    let memory_bytes = rss_before
        .zip(rss_after)
        .map(|(b, a)| a.saturating_sub(b))
        .filter(|&delta| delta >= 4096)
        .unwrap_or(theoretical);

    Ok(MetricsResult {
        fidelity,
        depth: stats.depth,
        gate_count: stats.gate_count,
        execution_time_ms: elapsed_ms,
        memory_bytes: Some(memory_bytes),
        backend_name: backend.name().to_string(),
    })
}
