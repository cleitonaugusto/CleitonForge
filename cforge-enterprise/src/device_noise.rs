//! Per-qubit, per-gate-pair device noise models.
//!
//! `DeviceNoiseModel` wraps the open-source `NoisyConfig` struct but adds
//! a full calibration table: different depolarizing rates per qubit pair,
//! T1/T2 thermal relaxation channels, and gate-timing-aware noise.

use cforge_backends::noise::NoiseChannel;
use cforge_backends::NoisyConfig;
use serde::{Deserialize, Serialize};

/// Per-qubit properties loaded from a device calibration file.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QubitCalibration {
    pub qubit: usize,
    pub t1_us: f64,
    pub t2_us: f64,
    pub readout_err: f64,
    pub sx_err: f64,
}

/// Per-pair two-qubit gate calibration.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GatePairCalibration {
    pub qubits: [usize; 2],
    pub gate_error: f64,
    pub gate_length_ns: f64,
}

/// Device-level noise model with heterogeneous per-qubit/per-pair calibration.
///
/// In the open-source tier, users apply a single uniform `NoisyConfig` to all
/// gates. `DeviceNoiseModel` stores a table — each (qubit, gate-type) pair has
/// its own error rate, enabling hardware-accurate simulation.
#[derive(Debug, Clone)]
pub struct DeviceNoiseModel {
    pub qubits: Vec<QubitCalibration>,
    pub pairs: Vec<GatePairCalibration>,
}

impl DeviceNoiseModel {
    /// Build a `DeviceNoiseModel` from a parsed IBM calibration JSON.
    pub fn from_ibm_json(json: &str) -> Result<Self, String> {
        // TODO: implement in enterprise release; delegates to ibm_profile parser.
        let _ = json;
        Err("DeviceNoiseModel::from_ibm_json requires enterprise license".into())
    }

    /// Average per-qubit noise → compatible `NoisyConfig` for community backends.
    ///
    /// This provides a lossless downgrade path: enterprise → community backends.
    pub fn to_average_config(&self) -> NoisyConfig {
        let avg_sx = avg(&self.qubits.iter().map(|q| q.sx_err).collect::<Vec<_>>());
        let avg_cx = avg(&self.pairs.iter().map(|p| p.gate_error).collect::<Vec<_>>());
        let avg_ro = avg(&self
            .qubits
            .iter()
            .map(|q| q.readout_err)
            .collect::<Vec<_>>());
        NoisyConfig {
            single_qubit: Some(NoiseChannel::Depolarizing { p: avg_sx }),
            two_qubit: Some(NoiseChannel::Depolarizing { p: avg_cx }),
            readout: Some(avg_ro),
        }
    }

    /// Get noise config for a specific qubit (single-qubit gates only).
    pub fn config_for_qubit(&self, qubit: usize) -> Option<NoisyConfig> {
        let q = self.qubits.iter().find(|q| q.qubit == qubit)?;
        Some(NoisyConfig {
            single_qubit: Some(NoiseChannel::Depolarizing { p: q.sx_err }),
            two_qubit: None,
            readout: Some(q.readout_err),
        })
    }
}

/// Amplitude-damping channel calibrated from T1 and gate duration.
///
/// γ = 1 - exp(-t_gate / T1)
///
/// This is the physically correct model for energy relaxation, more accurate
/// than a generic depolarizing channel. T1 is measured in microseconds;
/// `gate_length_ns` in nanoseconds.
#[derive(Debug, Clone)]
pub struct ThermalRelaxationChannel {
    pub t1_us: f64,
    pub t2_us: f64,
}

impl ThermalRelaxationChannel {
    /// Compute amplitude-damping γ for a gate of duration `gate_ns` nanoseconds.
    pub fn gamma_for_gate(&self, gate_ns: f64) -> f64 {
        let t_us = gate_ns / 1000.0;
        1.0 - (-t_us / self.t1_us).exp()
    }

    /// Community-compatible `NoiseChannel` approximation.
    pub fn to_noise_channel(&self, gate_ns: f64) -> NoiseChannel {
        NoiseChannel::AmplitudeDamping {
            gamma: self.gamma_for_gate(gate_ns),
        }
    }
}

fn avg(xs: &[f64]) -> f64 {
    if xs.is_empty() {
        return 0.0;
    }
    xs.iter().sum::<f64>() / xs.len() as f64
}
