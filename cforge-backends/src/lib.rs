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
//! - [`roqoqo_backend::RoqoqoBackend`] — uses roqoqo's gate matrices.
//! - [`q1tsim_backend::Q1tSimBackend`] — uses q1tsim's gate matrices.
//! - [`density_matrix::DensityMatrixBackend`] — exact noisy simulation via ρ.

pub mod trait_def;
pub mod statevector;
pub mod quantrs2_backend;
pub mod roqoqo_backend;
#[cfg(feature = "q1tsim")]
pub mod q1tsim_backend;
pub mod density_matrix;
pub mod noise;
pub mod noisy_backend;
mod sample;

pub use trait_def::{BackendError, DEFAULT_SEED, SimulationBackend, SimulationResult};
pub use statevector::NativeStateVectorBackend;
pub use quantrs2_backend::QuantRS2Backend;
pub use roqoqo_backend::RoqoqoBackend;
#[cfg(feature = "q1tsim")]
pub use q1tsim_backend::Q1tSimBackend;
pub use density_matrix::DensityMatrixBackend;
pub use noisy_backend::NoisyStatevectorBackend;
pub use noise::NoisyConfig;
