//! Pluggable simulation backends for CleitonForge.
//!
//! Every backend implements [`trait_def::SimulationBackend`] and
//! accepts a [`cforge_core::Circuit`] as input, so callers never need
//! to know which backend they're talking to.
//!
//! Available backends:
//! - [`statevector::NativeStateVectorBackend`] — custom state-vector
//!   simulation built entirely within CleitonForge using `num-complex`.
//! - [`quantrs2_backend::QuantRS2Backend`] — uses `quantrs2-core`'s
//!   gate matrix definitions as the source of truth for unitary matrices.

pub mod trait_def;
pub mod statevector;
pub mod quantrs2_backend;
mod sample;

pub use trait_def::{BackendError, SimulationBackend, SimulationResult};
pub use statevector::NativeStateVectorBackend;
pub use quantrs2_backend::QuantRS2Backend;
