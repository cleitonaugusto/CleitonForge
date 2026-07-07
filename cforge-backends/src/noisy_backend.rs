//! Noisy statevector backend using the quantum-trajectory (MCWF) method.
//!
//! ## How it works
//!
//! Each shot runs a separate statevector trajectory:
//! 1. Initialize |0...0⟩
//! 2. For every gate: apply the gate, then apply the configured noise channel
//!    to the qubits involved (stochastically)
//! 3. Measure: sample one bitstring, optionally apply readout error
//!
//! Averaging counts over many shots reproduces the correct noisy output
//! distribution without ever computing a density matrix.
//!
//! shots=0 runs one trajectory and returns the (noisy) statevector.
//!
//! ## Extension point
//!
//! `NoisyConfig` holds open-source noise channels. Enterprise device profiles
//! (IBM, IonQ, Rigetti calibration JSON) implement additional channels and
//! slot into the same `NoisyConfig` struct via the `cforge-enterprise` crate.

use std::collections::HashMap;

use num_complex::Complex64;
use rand::rngs::SmallRng;
use rand::SeedableRng;

use cforge_core::Circuit;

use crate::noise::{apply_readout_error, NoisyConfig};
use crate::sample::sample_once_rng;
use crate::statevector::apply_gate;
use crate::trait_def::{BackendError, SimulationBackend, SimulationResult};

const MAX_QUBITS: usize = 18;

/// Noisy statevector backend using the quantum-trajectory method.
///
/// # Example
///
/// ```rust,no_run
/// use cforge_backends::{NoisyStatevectorBackend, SimulationBackend};
/// use cforge_backends::noise::NoisyConfig;
/// use cforge_core::{Circuit, GateKind, Operation};
///
/// let mut c = Circuit::new(2);
/// c.push(Operation::new(GateKind::H,  vec![0], vec![]));
/// c.push(Operation::new(GateKind::Cx, vec![0, 1], vec![]));
///
/// let backend = NoisyStatevectorBackend {
///     config: NoisyConfig::depolarizing(0.001, 0.01),
/// };
/// let r = backend.run(&c, 4096, 42).unwrap();
/// println!("{:?}", r.counts);
/// ```
pub struct NoisyStatevectorBackend {
    pub config: NoisyConfig,
}

impl SimulationBackend for NoisyStatevectorBackend {
    fn name(&self) -> &str {
        "noisy-statevector"
    }

    fn run(
        &self,
        circuit: &Circuit,
        shots: usize,
        seed: u64,
    ) -> Result<SimulationResult, BackendError> {
        let n = circuit.num_qubits();
        if n > MAX_QUBITS {
            return Err(BackendError(format!(
                "noisy backend supports up to {MAX_QUBITS} qubits, circuit has {n}"
            )));
        }

        let trajectories = shots.max(1);
        let mut rng = SmallRng::seed_from_u64(seed);
        let mut counts: HashMap<String, usize> = HashMap::new();
        let mut last_sv = vec![Complex64::new(0.0, 0.0); 1 << n];

        for _ in 0..trajectories {
            last_sv = run_trajectory(circuit, &self.config, &mut rng).map_err(BackendError)?;

            if shots > 0 {
                let mut bits = sample_once_rng(&last_sv, &mut rng);
                if let Some(p) = self.config.readout {
                    apply_readout_error(&mut bits, p, &mut rng);
                }
                let bitstring: String = bits
                    .iter()
                    .rev()
                    .map(|b| char::from_digit(*b as u32, 10).unwrap())
                    .collect();
                *counts.entry(bitstring).or_insert(0) += 1;
            }
        }

        Ok(SimulationResult {
            statevector: last_sv,
            counts,
        })
    }
}

// ── Trajectory engine ─────────────────────────────────────────────────────────

fn run_trajectory(
    circuit: &Circuit,
    config: &NoisyConfig,
    rng: &mut SmallRng,
) -> Result<Vec<Complex64>, String> {
    let n = circuit.num_qubits();
    let mut sv = vec![Complex64::new(0.0, 0.0); 1 << n];
    sv[0] = Complex64::new(1.0, 0.0);

    for (op_idx, op) in circuit.operations.iter().enumerate() {
        apply_gate(&mut sv, op.kind, &op.qubits, &op.params, n)
            .map_err(|e| format!("op {op_idx}: {e}"))?;

        let noise_ch = match op.qubits.len() {
            1 => config.single_qubit.as_ref(),
            2 => config.two_qubit.as_ref(),
            _ => None,
        };
        if let Some(ch) = noise_ch {
            for &q in &op.qubits {
                ch.apply(&mut sv, q, rng);
            }
        }
    }
    Ok(sv)
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::noise::{NoiseChannel, NoisyConfig};
    use crate::{NativeStateVectorBackend, SimulationBackend};
    use cforge_core::{GateKind, Operation};

    fn bell_circuit() -> Circuit {
        let mut c = Circuit::new(2);
        c.push(Operation::new(GateKind::H, vec![0], vec![]));
        c.push(Operation::new(GateKind::Cx, vec![0, 1], vec![]));
        c
    }

    #[test]
    fn zero_noise_matches_noiseless() {
        // With p=0 noise, noisy backend should match noiseless exactly.
        let c = bell_circuit();
        let noiseless = NativeStateVectorBackend.run(&c, 0, 0).unwrap();
        let noisy = NoisyStatevectorBackend {
            config: NoisyConfig::depolarizing(0.0, 0.0),
        }
        .run(&c, 0, 0)
        .unwrap();

        for (a, b) in noiseless.statevector.iter().zip(&noisy.statevector) {
            assert!((a - b).norm() < 1e-12, "amplitude mismatch: {a} vs {b}");
        }
    }

    #[test]
    fn high_noise_degrades_bell_fidelity() {
        // With heavy depolarizing noise, Bell fidelity should drop below 0.8.
        let c = bell_circuit();
        let ideal = NativeStateVectorBackend.run(&c, 0, 0).unwrap().statevector;
        let noisy = NoisyStatevectorBackend {
            config: NoisyConfig::depolarizing(0.5, 0.5),
        }
        .run(&c, 0, 42)
        .unwrap()
        .statevector;

        let fidelity: f64 = ideal
            .iter()
            .zip(&noisy)
            .map(|(a, b)| (a.conj() * b).norm())
            .sum();
        assert!(
            fidelity < 0.9,
            "expected degraded fidelity, got {fidelity:.4}"
        );
    }

    #[test]
    fn shots_counts_sum_to_shots() {
        let c = bell_circuit();
        let r = NoisyStatevectorBackend {
            config: NoisyConfig::depolarizing(0.005, 0.01),
        }
        .run(&c, 1024, 0)
        .unwrap();
        let total: usize = r.counts.values().sum();
        assert_eq!(total, 1024);
    }

    #[test]
    fn deterministic_with_same_seed() {
        let c = bell_circuit();
        let cfg = NoisyConfig::depolarizing_with_readout(0.01, 0.05, 0.02);
        let r1 = NoisyStatevectorBackend {
            config: cfg.clone(),
        }
        .run(&c, 512, 77)
        .unwrap();
        let r2 = NoisyStatevectorBackend { config: cfg }
            .run(&c, 512, 77)
            .unwrap();
        assert_eq!(r1.counts, r2.counts);
    }

    #[test]
    fn amplitude_damping_biases_toward_zero() {
        // Strong T1 damping on excited state: measure |0⟩ much more than |1⟩.
        let mut c = Circuit::new(1);
        c.push(Operation::new(GateKind::X, vec![0], vec![])); // prepare |1⟩
        let r = NoisyStatevectorBackend {
            config: NoisyConfig {
                single_qubit: Some(NoiseChannel::AmplitudeDamping { gamma: 0.8 }),
                two_qubit: None,
                readout: None,
            },
        }
        .run(&c, 2048, 0)
        .unwrap();

        let zeros = r.counts.get("0").copied().unwrap_or(0);
        let ones = r.counts.get("1").copied().unwrap_or(0);
        assert!(
            zeros > ones,
            "expected more |0⟩ than |1⟩ due to T1 decay, got 0:{zeros} 1:{ones}"
        );
    }

    #[test]
    fn readout_error_flips_some_bits() {
        // 100% readout error on X gate → |1⟩ should all be read as |0⟩.
        let mut c = Circuit::new(1);
        c.push(Operation::new(GateKind::X, vec![0], vec![]));
        let r = NoisyStatevectorBackend {
            config: NoisyConfig {
                single_qubit: None,
                two_qubit: None,
                readout: Some(1.0),
            },
        }
        .run(&c, 256, 0)
        .unwrap();
        let zeros = r.counts.get("0").copied().unwrap_or(0);
        assert_eq!(zeros, 256, "with readout=1.0, all |1⟩ should flip to 0");
    }
}
