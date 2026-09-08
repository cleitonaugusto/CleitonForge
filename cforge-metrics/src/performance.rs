//! Wall-clock timing and peak-memory measurement around backend calls.

use cforge_backends::{BackendError, SimulationBackend, SimulationResult};
use cforge_core::{Circuit, MetricsResult};
use std::time::Instant;

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
        counts: result.counts,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use cforge_backends::{NativeStateVectorBackend, DEFAULT_SEED};
    use cforge_core::{GateKind, Operation};

    fn bell() -> Circuit {
        let mut c = Circuit::new(2);
        c.push(Operation::new(GateKind::H, vec![0], vec![]));
        c.push(Operation::new(GateKind::Cx, vec![0, 1], vec![]));
        c
    }

    /// The backend produces counts and `measure` used to drop them, so the CLI
    /// could run shots and never say what was measured.
    #[test]
    fn shots_produce_counts_that_sum_to_the_shot_count() {
        let m = measure(&NativeStateVectorBackend, &bell(), 2048, DEFAULT_SEED, None).unwrap();
        assert_eq!(m.counts.values().sum::<usize>(), 2048);
        // A Bell state only ever lands on |00> or |11>.
        let mut seen: Vec<&str> = m.counts.keys().map(String::as_str).collect();
        seen.sort_unstable();
        assert_eq!(seen, ["00", "11"]);
    }

    /// Without shots there is nothing to count, and an empty map says that
    /// more honestly than zeros would.
    #[test]
    fn no_shots_means_no_counts() {
        let m = measure(&NativeStateVectorBackend, &bell(), 0, DEFAULT_SEED, None).unwrap();
        assert!(m.counts.is_empty());
    }

    /// The seed is documented as making a run reproducible. That is a claim
    /// about behaviour, so it gets a test.
    #[test]
    fn the_same_seed_gives_the_same_counts() {
        let a = measure(&NativeStateVectorBackend, &bell(), 2048, 42, None).unwrap();
        let b = measure(&NativeStateVectorBackend, &bell(), 2048, 42, None).unwrap();
        assert_eq!(a.counts, b.counts);
    }
}
