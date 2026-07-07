//! `cforge-enterprise` — device-specific noise calibration layer.
//!
//! ## Open-core boundary
//!
//! Community (`cforge-backends`):
//! - `NoiseChannel` trait + 4 noise models (Depolarizing, BitFlip, PhaseFlip,
//!   AmplitudeDamping)
//! - `NoisyConfig` builder (uniform noise across all gates)
//! - `NoisyStatevectorBackend` (Monte Carlo trajectories)
//! - `DensityMatrixBackend` (exact noisy simulation)
//!
//! Enterprise (this crate, private/commercial):
//! - `DeviceNoiseModel` — per-qubit, per-gate-pair calibration tables
//! - `LiveCalibration` — live API pull from IBM/IonQ/Google
//! - `ZneEstimator` — Zero-Noise Extrapolation for error mitigation
//! - `ReadoutMitigationMatrix` — measurement error correction via M-inverse
//! - `ThermalRelaxationChannel` — T1/T2-accurate depolarizing from timing
//!
//! ## Why private
//!
//! Device calibration data requires API keys and commercial partnerships.
//! Mitigation algorithms are differentiating IP. The open-source layer
//! (`NoiseChannel` trait) ensures enterprise extensions are drop-in replacements
//! for community noise models — users can switch without changing their circuit code.

pub mod device_noise;
pub mod readout_mitigation;
pub mod zne;

pub use device_noise::{DeviceNoiseModel, ThermalRelaxationChannel};
pub use readout_mitigation::ReadoutMitigationMatrix;
pub use zne::ZneEstimator;
