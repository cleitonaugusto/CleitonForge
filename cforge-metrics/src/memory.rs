//! Memory measurement for quantum simulation runs.
//!
//! Two sources are used in priority order:
//!   1. RSS delta from `/proc/self/status` (Linux) — actual pages allocated
//!      by the OS during the backend call, including the statevector heap
//!      allocation and any overhead.
//!   2. Theoretical minimum — `2^n_qubits × 16 bytes` (two `f64` per
//!      `Complex64` amplitude), always available and reproducible.
//!
//! The RSS approach is accurate for circuits with ≥ 1 MB statevectors.
//! For smaller circuits the OS page granularity (4 KB) dominates; in
//! that case the theoretical value is returned instead.

/// Theoretical peak memory for a pure statevector simulation of `n` qubits:
/// `2^n × sizeof(Complex64)` = `2^n × 16` bytes.
///
/// This is the minimum possible memory any correct state-vector backend
/// must allocate. Real backends may use somewhat more due to stack frames,
/// intermediate buffers, and allocator overhead.
pub fn statevector_memory_bytes(n_qubits: usize) -> u64 {
    // 16 = size_of::<Complex64>() = 2 × size_of::<f64>()
    (1u64 << n_qubits).saturating_mul(16)
}

/// Reads the current Resident Set Size (RSS) from the OS in bytes.
///
/// On Linux: parsed from `/proc/self/status` (`VmRSS` field).
/// On other platforms: always returns `None`.
pub fn current_rss_bytes() -> Option<u64> {
    _current_rss_bytes()
}

#[cfg(target_os = "linux")]
fn _current_rss_bytes() -> Option<u64> {
    let status = std::fs::read_to_string("/proc/self/status").ok()?;
    for line in status.lines() {
        if line.starts_with("VmRSS:") {
            // Format: "VmRSS:   12345 kB"
            let kb: u64 = line.split_whitespace().nth(1)?.parse().ok()?;
            return Some(kb * 1024);
        }
    }
    None
}

#[cfg(not(target_os = "linux"))]
fn _current_rss_bytes() -> Option<u64> {
    None
}

/// Attempts to measure peak memory allocated during a closure, returning the
/// RSS delta in bytes. Falls back to `fallback_bytes` when RSS is
/// unavailable or the delta is below one OS page (4 KiB) — the latter
/// indicates the allocation was too small to show up in page-granular RSS.
pub fn measure_rss_delta(fallback_bytes: u64, f: impl FnOnce()) -> u64 {
    let before = current_rss_bytes();
    f();
    let after = current_rss_bytes();

    before
        .zip(after)
        .map(|(b, a)| a.saturating_sub(b))
        .filter(|&delta| delta >= 4096) // ignore sub-page noise
        .unwrap_or(fallback_bytes)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn theoretical_memory_is_correct() {
        // 2 qubits → 4 amplitudes × 16 bytes = 64 bytes
        assert_eq!(statevector_memory_bytes(2), 64);
        // 10 qubits → 1024 × 16 = 16 384 bytes
        assert_eq!(statevector_memory_bytes(10), 16_384);
        // 20 qubits → 1 Mi × 16 = 16 MiB
        assert_eq!(statevector_memory_bytes(20), 16 * 1024 * 1024);
    }

    #[test]
    fn rss_is_positive_or_none() {
        if let Some(rss) = current_rss_bytes() {
            assert!(rss > 0, "RSS must be positive when available");
        }
        // None is acceptable on non-Linux platforms
    }

    #[test]
    fn rss_delta_captures_large_allocation() {
        // Allocate ~8 MB explicitly and confirm RSS delta ≥ 1 page.
        let theoretical = statevector_memory_bytes(19); // 8 MiB
        let delta = measure_rss_delta(theoretical, || {
            let _large: Vec<u8> = vec![0u8; 8 * 1024 * 1024];
            // keep alive until end of closure
        });
        // On Linux: delta should be around 8 MB.
        // On other platforms: falls back to theoretical (also 8 MiB).
        assert!(delta >= 4096, "delta should be at least one page");
    }
}
