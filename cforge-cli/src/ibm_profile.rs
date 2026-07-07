//! Parse IBM Quantum JSON calibration files → `NoisyConfig`.
//!
//! IBM calibration JSON is public (available from IBM Quantum Experience).
//! This parser converts the standard format into a `NoisyConfig` that can
//! be passed to `NoisyStatevectorBackend` or `DensityMatrixBackend`.
//!
//! ## Example
//!
//! ```ignore
//! let json = std::fs::read_to_string("examples/ibm_noise_profile.json")?;
//! let cfg = parse_ibm_profile(&json)?;
//! let backend = NoisyStatevectorBackend { config: cfg };
//! ```

use cforge_backends::NoisyConfig;

#[derive(Debug)]
pub struct IbmProfile {
    pub backend: String,
    pub avg_sx_error: f64,
    pub avg_cx_error: f64,
    pub avg_readout_error: f64,
}

impl IbmProfile {
    pub fn to_noisy_config(&self) -> NoisyConfig {
        NoisyConfig::depolarizing_with_readout(
            self.avg_sx_error,
            self.avg_cx_error,
            self.avg_readout_error,
        )
    }
}

/// Parse a subset of IBM calibration JSON into an `IbmProfile`.
///
/// The parser is deliberately simple — no serde dependency required.
/// It reads average SX, CX, and readout error rates from the JSON.
pub fn parse_ibm_profile(json: &str) -> Result<IbmProfile, String> {
    // Use serde_json (already a dep in cforge-cli) for proper parsing
    let v: serde_json::Value =
        serde_json::from_str(json).map_err(|e| format!("JSON parse error: {e}"))?;

    let backend = v["backend"].as_str().unwrap_or("unknown").to_string();

    // Average SX error across all qubits
    let avg_sx_error = average_gate_errors(&v["gate_errors"]["sx"]);
    // Average CX error across all pairs
    let avg_cx_error = average_gate_errors(&v["gate_errors"]["cx"]);
    // Average readout error across all qubits
    let avg_readout_error = {
        let arr = v["qubit_props"].as_array().ok_or("missing qubit_props")?;
        let sum: f64 = arr.iter().filter_map(|q| q["readout_err"].as_f64()).sum();
        sum / arr.len().max(1) as f64
    };

    if avg_sx_error == 0.0 && avg_cx_error == 0.0 {
        return Err("no gate error data found in profile".to_string());
    }

    Ok(IbmProfile {
        backend,
        avg_sx_error,
        avg_cx_error,
        avg_readout_error,
    })
}

fn average_gate_errors(arr: &serde_json::Value) -> f64 {
    let entries = match arr.as_array() {
        Some(v) => v,
        None => return 0.0,
    };
    let (sum, count) = entries
        .iter()
        .filter_map(|e| e["gate_error"].as_f64())
        .fold((0.0, 0usize), |(s, n), p| (s + p, n + 1));
    if count == 0 {
        0.0
    } else {
        sum / count as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_nairobi_snapshot() {
        let json = include_str!("../../examples/ibm_noise_profile.json");
        let p = parse_ibm_profile(json).expect("parse failed");
        assert_eq!(p.backend, "ibm_nairobi");
        assert!(
            p.avg_sx_error > 1e-5 && p.avg_sx_error < 0.01,
            "sx error = {}",
            p.avg_sx_error
        );
        assert!(
            p.avg_cx_error > 0.001 && p.avg_cx_error < 0.05,
            "cx error = {}",
            p.avg_cx_error
        );
        assert!(
            p.avg_readout_error > 1e-4 && p.avg_readout_error < 0.1,
            "readout = {}",
            p.avg_readout_error
        );
    }

    #[test]
    fn profile_to_noisy_config_valid() {
        let json = include_str!("../../examples/ibm_noise_profile.json");
        let p = parse_ibm_profile(json).unwrap();
        let cfg = p.to_noisy_config();
        assert!(cfg.single_qubit.is_some());
        assert!(cfg.two_qubit.is_some());
        assert!(cfg.readout.is_some());
    }
}
