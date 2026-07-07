//! NISQ noise channels for stochastic statevector simulation.
//!
//! Each channel implements [`NoiseChannel::apply`], which mutates a statevector
//! in-place using the quantum-trajectory (Monte Carlo wavefunction) method.
//! The noisy backend calls this after every gate; averaging over many shots
//! reproduces the correct noisy probability distribution.
//!
//! ## Channels available
//!
//! | Channel | Description |
//! |---|---|
//! | `Depolarizing` | I w.p. (1-p), X/Y/Z each w.p. p/3 |
//! | `BitFlip` | X w.p. p |
//! | `PhaseFlip` | Z w.p. p |
//! | `AmplitudeDamping` | Quantum jump: energy decay \|1⟩→\|0⟩ w.p. γ·P(\|1⟩) |
//!
//! ## SaaS extension point
//!
//! Production noise profiles (device-specific gate error rates, T1/T2 from
//! calibration JSON, correlated 2-qubit errors) are implemented as additional
//! `NoiseChannel` variants in the private `cforge-enterprise` crate that
//! depends on this public API.

use num_complex::Complex64;
use rand::{Rng, RngExt};

// ── Internal statevector primitive ────────────────────────────────────────────

type U2 = [[Complex64; 2]; 2];

fn apply1_sv(sv: &mut [Complex64], k: usize, u: &U2) {
    let stride = 1 << k;
    let len = sv.len();
    let mut i = 0;
    while i < len {
        for j in i..(i + stride) {
            let a = sv[j];
            let b = sv[j + stride];
            sv[j] = u[0][0] * a + u[0][1] * b;
            sv[j + stride] = u[1][0] * a + u[1][1] * b;
        }
        i += 2 * stride;
    }
}

fn renorm(sv: &mut [Complex64]) {
    let norm: f64 = sv.iter().map(|a| a.norm_sqr()).sum::<f64>().sqrt();
    if norm > 1e-14 {
        for a in sv.iter_mut() {
            *a /= norm;
        }
    }
}

const X: U2 = [
    [Complex64::new(0.0, 0.0), Complex64::new(1.0, 0.0)],
    [Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)],
];
const Y: U2 = [
    [Complex64::new(0.0, 0.0), Complex64::new(0.0, -1.0)],
    [Complex64::new(0.0, 1.0), Complex64::new(0.0, 0.0)],
];
const Z: U2 = [
    [Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)],
    [Complex64::new(0.0, 0.0), Complex64::new(-1.0, 0.0)],
];

// ── Public channel type ───────────────────────────────────────────────────────

/// A noise channel applied after each gate in the quantum-trajectory method.
///
/// Construct via the named constructors. Configure a [`NoisyConfig`] and pass
/// it to [`NoisyStatevectorBackend`][crate::noisy_backend::NoisyStatevectorBackend].
#[derive(Clone, Debug)]
pub enum NoiseChannel {
    /// Symmetric depolarizing: applies X, Y, or Z each with probability p/3.
    Depolarizing { p: f64 },
    /// Bit-flip: applies X with probability p.
    BitFlip { p: f64 },
    /// Phase-flip: applies Z with probability p.
    PhaseFlip { p: f64 },
    /// Amplitude damping (energy decay |1⟩→|0⟩).
    ///
    /// `gamma` is the single-step decay probability. The jump probability for a
    /// given trajectory step is `gamma * P(qubit=|1⟩)`.
    AmplitudeDamping { gamma: f64 },
}

impl NoiseChannel {
    /// Apply this channel stochastically to `qubit` in the statevector `sv`.
    pub fn apply(&self, sv: &mut [Complex64], qubit: usize, rng: &mut impl Rng) {
        match *self {
            NoiseChannel::Depolarizing { p } => apply_depolarizing(sv, qubit, p, rng),
            NoiseChannel::BitFlip { p } => {
                if rng.random::<f64>() < p {
                    apply1_sv(sv, qubit, &X);
                }
            }
            NoiseChannel::PhaseFlip { p } => {
                if rng.random::<f64>() < p {
                    apply1_sv(sv, qubit, &Z);
                }
            }
            NoiseChannel::AmplitudeDamping { gamma } => {
                apply_amplitude_damping(sv, qubit, gamma, rng);
            }
        }
    }
}

// ── Channel implementations ───────────────────────────────────────────────────

fn apply_depolarizing(sv: &mut [Complex64], k: usize, p: f64, rng: &mut impl Rng) {
    let r = rng.random::<f64>();
    if r < p / 3.0 {
        apply1_sv(sv, k, &X);
    } else if r < 2.0 * p / 3.0 {
        apply1_sv(sv, k, &Y);
    } else if r < p {
        apply1_sv(sv, k, &Z);
    }
    // r >= p → no error (probability 1-p)
}

/// Quantum-jump amplitude damping via the MCWF method.
///
/// Kraus operators:
///   K0 = [[1, 0], [0, √(1-γ)]]   (no decay)
///   K1 = [[0, √γ], [0, 0]]        (decay: |1⟩→|0⟩)
fn apply_amplitude_damping(sv: &mut [Complex64], k: usize, gamma: f64, rng: &mut impl Rng) {
    let bit = 1 << k;

    // P(|1⟩ for qubit k) — probability of jump
    let p1: f64 = sv
        .iter()
        .enumerate()
        .filter(|(i, _)| i & bit != 0)
        .map(|(_, a)| a.norm_sqr())
        .sum();

    let p_jump = gamma * p1;

    if rng.random::<f64>() < p_jump {
        // Apply K1: |1⟩ component collapses to |0⟩
        for i in 0..sv.len() {
            if i & bit == 0 {
                sv[i] = sv[i | bit] * gamma.sqrt();
                sv[i | bit] = Complex64::new(0.0, 0.0);
            }
        }
    } else {
        // Apply K0: damp |1⟩ amplitude
        let s = (1.0 - gamma).sqrt();
        for (i, amp) in sv.iter_mut().enumerate() {
            if i & bit != 0 {
                *amp *= s;
            }
        }
    }
    renorm(sv);
}

// ── Noise configuration ───────────────────────────────────────────────────────

/// Noise configuration for [`NoisyStatevectorBackend`][crate::noisy_backend::NoisyStatevectorBackend].
///
/// # Common NISQ device parameters (order-of-magnitude)
///
/// | Parameter | Typical value |
/// |---|---|
/// | Single-qubit gate depolarizing | p = 0.001 (0.1%) |
/// | Two-qubit gate depolarizing | p = 0.01 (1%) |
/// | Readout error | p = 0.01–0.05 |
/// | T1 amplitude damping per gate | gamma ≈ gate_time / T1 |
///
/// For certified device-specific profiles (IBM, IonQ, Rigetti, native QPU
/// calibration JSON), see the `cforge-enterprise` crate.
#[derive(Clone, Debug, Default)]
pub struct NoisyConfig {
    /// Noise applied after each single-qubit gate.
    pub single_qubit: Option<NoiseChannel>,
    /// Noise applied to each qubit after each two-qubit gate.
    pub two_qubit: Option<NoiseChannel>,
    /// Classical bit-flip probability at readout (symmetric).
    pub readout: Option<f64>,
}

impl NoisyConfig {
    /// Uniform depolarizing noise: `p1` after 1q gates, `p2` after 2q gates.
    pub fn depolarizing(p1: f64, p2: f64) -> Self {
        Self {
            single_qubit: Some(NoiseChannel::Depolarizing { p: p1 }),
            two_qubit: Some(NoiseChannel::Depolarizing { p: p2 }),
            readout: None,
        }
    }

    /// Depolarizing + classical readout error.
    pub fn depolarizing_with_readout(p1: f64, p2: f64, readout: f64) -> Self {
        Self {
            single_qubit: Some(NoiseChannel::Depolarizing { p: p1 }),
            two_qubit: Some(NoiseChannel::Depolarizing { p: p2 }),
            readout: Some(readout),
        }
    }

    /// Amplitude damping only (T1 decay per gate step).
    pub fn amplitude_damping(gamma: f64) -> Self {
        Self {
            single_qubit: Some(NoiseChannel::AmplitudeDamping { gamma }),
            two_qubit: Some(NoiseChannel::AmplitudeDamping { gamma }),
            readout: None,
        }
    }
}

// ── Readout noise ─────────────────────────────────────────────────────────────

/// Applies symmetric bit-flip readout error to a sampled bitstring.
pub(crate) fn apply_readout_error(bits: &mut [u8], p: f64, rng: &mut impl Rng) {
    for b in bits.iter_mut() {
        if rng.random::<f64>() < p {
            *b ^= 1;
        }
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::SmallRng;
    use rand::SeedableRng;

    fn bell_sv() -> Vec<Complex64> {
        let s = 1.0 / std::f64::consts::SQRT_2;
        // |Φ+⟩ = (|00⟩ + |11⟩)/√2 on qubits 1,0
        let mut sv = vec![Complex64::new(0.0, 0.0); 4];
        sv[0b00] = Complex64::new(s, 0.0);
        sv[0b11] = Complex64::new(s, 0.0);
        sv
    }

    #[test]
    fn zero_noise_identity() {
        let mut sv = bell_sv();
        let orig = sv.clone();
        let mut rng = SmallRng::seed_from_u64(0);
        // p=0 → never applies an error
        NoiseChannel::Depolarizing { p: 0.0 }.apply(&mut sv, 0, &mut rng);
        for (a, b) in sv.iter().zip(&orig) {
            assert!((a - b).norm() < 1e-14);
        }
    }

    #[test]
    fn depolarizing_full_noise_always_errors() {
        // With p=1 on |0⟩: X flips to |1⟩, Y flips to i|1⟩, Z is a global phase (invisible).
        // Expected visible changes: ~2/3 of shots (X and Y each with prob 1/3).
        let orig = vec![Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)];
        let mut changed = 0;
        let mut rng = SmallRng::seed_from_u64(42);
        for _ in 0..1200 {
            let mut s = orig.clone();
            NoiseChannel::Depolarizing { p: 1.0 }.apply(&mut s, 0, &mut rng);
            // Check if qubit flipped (|1⟩ amplitude increased)
            if s[1].norm() > 0.5 {
                changed += 1;
            }
        }
        // ~2/3 of shots should flip (X or Y applied); allow ±10%
        assert!(
            changed > 650 && changed < 950,
            "expected ~800 flips/1200 shots, got {changed}"
        );
    }

    #[test]
    fn amplitude_damping_decays_excited_state() {
        // |1⟩ state: with γ=1 always decays to |0⟩
        let mut sv = vec![Complex64::new(0.0, 0.0), Complex64::new(1.0, 0.0)];
        let mut rng = SmallRng::seed_from_u64(0);
        NoiseChannel::AmplitudeDamping { gamma: 1.0 }.apply(&mut sv, 0, &mut rng);
        assert!(sv[1].norm() < 1e-10, "|1⟩ should have decayed");
        assert!((sv[0].norm() - 1.0).abs() < 1e-10, "|0⟩ should be 1");
    }

    #[test]
    fn amplitude_damping_preserves_ground_state() {
        // |0⟩ never decays (no energy to lose)
        let mut sv = vec![Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)];
        let orig = sv.clone();
        let mut rng = SmallRng::seed_from_u64(0);
        for _ in 0..100 {
            let mut s = orig.clone();
            NoiseChannel::AmplitudeDamping { gamma: 0.5 }.apply(&mut s, 0, &mut rng);
            assert!((s[0].norm() - 1.0).abs() < 1e-10);
        }
    }

    #[test]
    fn readout_error_flips_bits() {
        let mut bits = vec![0u8, 0, 0, 0];
        let mut rng = SmallRng::seed_from_u64(99);
        apply_readout_error(&mut bits, 1.0, &mut rng);
        assert!(
            bits.iter().all(|&b| b == 1),
            "all bits should flip with p=1"
        );
    }
}
