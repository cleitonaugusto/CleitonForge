//! Standardized, backend-agnostic metrics computation for CleitonForge.
//!
//! All metrics operate on the canonical [`cforge_core::Circuit`] IR —
//! never on backend-native representations — so the numbers are
//! comparable regardless of which backend produced them.

pub mod circuit_stats;
pub mod fidelity;
pub mod performance;

pub use circuit_stats::{compute_stats, CircuitStats};
pub use fidelity::statevector_fidelity;
pub use performance::measure;
