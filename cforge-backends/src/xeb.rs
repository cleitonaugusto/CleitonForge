//! Linear cross-entropy benchmarking (XEB) and conjugated reference
//! backends for benchmark-blindness experiments.
//!
//! XEB is the fidelity proxy used in Google's 2019 quantum-supremacy
//! experiment (Arute et al., Nature 574): sample bitstrings from the
//! device, look up the ideal probability of each sample, and compute
//!
//! ```text
//! F_XEB = 2^n · E_x[ p_ideal(x) ] − 1
//! ```
//!
//! Because F_XEB depends only on the ideal *probabilities* |⟨x|C|0⟩|²,
//! it inherits the conjugation blindness proved in the QGCS paper: a
//! simulator that implements every gate as its complex conjugate U*
//! produces exactly the same XEB score as a correct one. The backends in
//! this module make that failure mode concrete and testable.

use std::collections::HashMap;

use num_complex::Complex64;

use cforge_core::{Circuit, GateKind, Operation};

use crate::sample::{sample_counts, Lcg64};
use crate::statevector::apply_gate_ext;
use crate::trait_def::{BackendError, SimulationBackend, SimulationResult};

const MAX_QUBITS: usize = 22;

/// Which gates the conjugation bug affects.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConjugationScope {
    /// Every gate is implemented as its complex conjugate — a consistent
    /// library-wide convention error (e.g. a global i → −i sign flip).
    /// This is a valid antiunitary image of quantum mechanics, so *no*
    /// benchmark that samples a circuit from |0…0⟩ can detect it.
    AllGates,
    /// Only the diagonal rotation family (Rz, Crz, Phase, Cp) has its
    /// rotation sign inverted — the class of bug found in quantrs2
    /// v0.2.0. On circuits whose remaining gates are all real (QV's
    /// Rz/Ry/CX basis, textbook QPE's H/X/Swap), this is equivalent to
    /// `AllGates` and equally invisible; it becomes detectable only when
    /// combined with a correctly-implemented complex gate such as S.
    RotationSignOnly,
}

/// A statevector simulator with a deliberate, physically-consistent
/// conjugation bug, used as the "wrong but plausible" device in
/// blindness experiments.
pub struct ConjugatedStateVectorBackend {
    pub scope: ConjugationScope,
}

impl ConjugatedStateVectorBackend {
    pub fn new(scope: ConjugationScope) -> Self {
        Self { scope }
    }

    fn conjugates(&self, kind: GateKind) -> bool {
        match self.scope {
            ConjugationScope::AllGates => true,
            ConjugationScope::RotationSignOnly => matches!(
                kind,
                GateKind::Rz | GateKind::Crz | GateKind::Phase | GateKind::Cp
            ),
        }
    }
}

impl SimulationBackend for ConjugatedStateVectorBackend {
    fn name(&self) -> &str {
        match self.scope {
            ConjugationScope::AllGates => "statevector-conjugated",
            ConjugationScope::RotationSignOnly => "statevector-rz-sign-bug",
        }
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
                "conjugated backend supports up to {MAX_QUBITS} qubits, circuit has {n}"
            )));
        }

        let mut sv = vec![Complex64::new(0.0, 0.0); 1 << n];
        sv[0] = Complex64::new(1.0, 0.0);

        for (op_idx, op) in circuit.operations.iter().enumerate() {
            apply_gate_ext(
                &mut sv,
                op.kind,
                &op.qubits,
                &op.params,
                n,
                self.conjugates(op.kind),
            )
            .map_err(|e| BackendError(format!("operation {op_idx}: {e}")))?;
        }

        let counts = if shots > 0 {
            sample_counts(&sv, shots, seed)
        } else {
            HashMap::new()
        };

        Ok(SimulationResult {
            statevector: sv,
            counts,
        })
    }
}

// ── XEB scores ───────────────────────────────────────────────────────────────

/// Exact (infinite-shot) linear XEB score:
/// `F = 2^n · Σ_x p_ideal(x)·p_device(x) − 1`.
///
/// For a perfect device on a Porter-Thomas-distributed random circuit
/// this approaches 1; for a device sampling uniformly it is 0.
pub fn xeb_score_exact(
    circuit: &Circuit,
    ideal: &dyn SimulationBackend,
    device: &dyn SimulationBackend,
) -> Result<f64, BackendError> {
    let p_ideal = ideal.run(circuit, 0, 0)?.probabilities();
    let p_device = device.run(circuit, 0, 0)?.probabilities();
    if p_ideal.len() != p_device.len() {
        return Err(BackendError(format!(
            "backend dimension mismatch: {} vs {}",
            p_ideal.len(),
            p_device.len()
        )));
    }
    let dim = p_ideal.len() as f64;
    let overlap: f64 = p_ideal.iter().zip(&p_device).map(|(a, b)| a * b).sum();
    Ok(dim * overlap - 1.0)
}

/// Sampled linear XEB score, as measured in a real experiment:
/// draw `shots` bitstrings from `device`, average the *ideal*
/// probability of each sample, then `F = 2^n · mean − 1`.
pub fn xeb_score_sampled(
    circuit: &Circuit,
    ideal: &dyn SimulationBackend,
    device: &dyn SimulationBackend,
    shots: usize,
    seed: u64,
) -> Result<f64, BackendError> {
    let p_ideal = ideal.run(circuit, 0, 0)?.probabilities();
    let counts = device.run(circuit, shots, seed)?.counts;
    let total: usize = counts.values().sum();
    if total == 0 {
        return Err(BackendError("device returned no samples".into()));
    }
    let mut acc = 0.0;
    for (bits, &count) in &counts {
        let idx = usize::from_str_radix(bits, 2)
            .map_err(|e| BackendError(format!("bad bitstring {bits:?}: {e}")))?;
        let p = p_ideal
            .get(idx)
            .ok_or_else(|| BackendError(format!("bitstring {bits:?} out of range")))?;
        acc += p * count as f64;
    }
    let dim = p_ideal.len() as f64;
    Ok(dim * (acc / total as f64) - 1.0)
}

/// Heavy-output-generation probability, the Quantum Volume figure of
/// merit: the fraction of device samples whose ideal probability
/// exceeds the median ideal probability. Ideal ≈ 0.85 (asymptotically
/// (1+ln 2)/2) for random circuits; 0.5 for uniform sampling.
pub fn heavy_output_probability(
    circuit: &Circuit,
    ideal: &dyn SimulationBackend,
    device: &dyn SimulationBackend,
    shots: usize,
    seed: u64,
) -> Result<f64, BackendError> {
    let p_ideal = ideal.run(circuit, 0, 0)?.probabilities();
    let mut sorted = p_ideal.clone();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let median = if sorted.len() % 2 == 0 {
        (sorted[sorted.len() / 2 - 1] + sorted[sorted.len() / 2]) / 2.0
    } else {
        sorted[sorted.len() / 2]
    };

    let counts = device.run(circuit, shots, seed)?.counts;
    let total: usize = counts.values().sum();
    if total == 0 {
        return Err(BackendError("device returned no samples".into()));
    }
    let mut heavy = 0usize;
    for (bits, &count) in &counts {
        let idx = usize::from_str_radix(bits, 2)
            .map_err(|e| BackendError(format!("bad bitstring {bits:?}: {e}")))?;
        if p_ideal[idx] > median {
            heavy += count;
        }
    }
    Ok(heavy as f64 / total as f64)
}

// ── Random circuit generation ────────────────────────────────────────────────

/// Generates a Quantum-Volume-style random circuit over the Rz/Ry/CX
/// basis: `depth` layers, each applying a random Euler rotation
/// (Rz·Ry·Rz) to every qubit followed by a chain of CX entanglers with
/// a randomly rotated starting offset. Deterministic in `seed`.
///
/// Every non-Rz gate in this basis is real, which is exactly why a
/// rotation-sign bug behaves like a global conjugation on these
/// circuits and slips past QV and XEB.
pub fn random_qv_circuit(num_qubits: usize, depth: usize, seed: u64) -> Circuit {
    let mut rng = Lcg64::new(seed);
    let mut angle = |scale: f64| (rng.next_f64() * 2.0 - 1.0) * scale;

    let mut c = Circuit::new(num_qubits);
    for layer in 0..depth {
        for q in 0..num_qubits {
            c.push(Operation::new(
                GateKind::Rz,
                vec![q],
                vec![angle(std::f64::consts::PI)],
            ));
            c.push(Operation::new(
                GateKind::Ry,
                vec![q],
                vec![angle(std::f64::consts::PI)],
            ));
            c.push(Operation::new(
                GateKind::Rz,
                vec![q],
                vec![angle(std::f64::consts::PI)],
            ));
        }
        if num_qubits > 1 {
            for q in 0..num_qubits - 1 {
                let (a, b) = ((q + layer) % num_qubits, (q + layer + 1) % num_qubits);
                if a != b {
                    c.push(Operation::new(GateKind::Cx, vec![a, b], vec![]));
                }
            }
        }
    }
    c
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::statevector::NativeStateVectorBackend;

    const SHOTS: usize = 20_000;

    /// The conjugated backend must produce the element-wise complex
    /// conjugate of the native statevector — the exact premise of the
    /// blindness theorem.
    #[test]
    fn conjugated_backend_conjugates_the_state() {
        let c = random_qv_circuit(4, 4, 7);
        let native = NativeStateVectorBackend.run(&c, 0, 0).unwrap();
        let conj = ConjugatedStateVectorBackend::new(ConjugationScope::AllGates)
            .run(&c, 0, 0)
            .unwrap();
        for (a, b) in native.statevector.iter().zip(&conj.statevector) {
            assert!((a.conj() - b).norm() < 1e-12);
        }
    }

    /// Corollary (XEB blindness): the conjugated device scores exactly
    /// what the correct device scores.
    #[test]
    fn xeb_is_blind_to_conjugation() {
        for seed in 0..5 {
            let c = random_qv_circuit(5, 5, seed);
            let ideal = NativeStateVectorBackend;
            let correct = NativeStateVectorBackend;
            let wrong = ConjugatedStateVectorBackend::new(ConjugationScope::AllGates);

            let f_correct = xeb_score_exact(&c, &ideal, &correct).unwrap();
            let f_wrong = xeb_score_exact(&c, &ideal, &wrong).unwrap();
            assert!(
                (f_correct - f_wrong).abs() < 1e-12,
                "seed {seed}: {f_correct} vs {f_wrong}"
            );
            // And the score is meaningful (near 1 for random circuits),
            // not trivially zero.
            assert!(f_correct > 0.3, "seed {seed}: degenerate score {f_correct}");
        }
    }

    /// The Rz-only bug (quantrs2 class) is equally invisible on QV
    /// circuits because every other gate in the basis is real.
    #[test]
    fn rotation_sign_bug_is_blind_on_qv_basis() {
        let c = random_qv_circuit(5, 5, 42);
        let ideal = NativeStateVectorBackend;
        let wrong = ConjugatedStateVectorBackend::new(ConjugationScope::RotationSignOnly);
        let f_ideal = xeb_score_exact(&c, &ideal, &ideal).unwrap();
        let f_wrong = xeb_score_exact(&c, &ideal, &wrong).unwrap();
        assert!((f_ideal - f_wrong).abs() < 1e-12);
    }

    /// Minimal probability-level witness for the Rz-only bug: H·Rz·S·H
    /// mixes the bugged rotation with a *correct* complex gate (S), so
    /// the output probabilities differ — this bug is detectable in
    /// principle, just not by benchmarks over real-plus-Rz bases.
    #[test]
    fn rotation_sign_bug_is_visible_to_a_crafted_witness() {
        let theta = 1.0;
        let mut c = Circuit::new(1);
        c.push(Operation::new(GateKind::H, vec![0], vec![]));
        c.push(Operation::new(GateKind::Rz, vec![0], vec![theta]));
        c.push(Operation::new(GateKind::S, vec![0], vec![]));
        c.push(Operation::new(GateKind::H, vec![0], vec![]));

        let good = NativeStateVectorBackend.run(&c, 0, 0).unwrap();
        let bad = ConjugatedStateVectorBackend::new(ConjugationScope::RotationSignOnly)
            .run(&c, 0, 0)
            .unwrap();
        let dp = (good.probabilities()[0] - bad.probabilities()[0]).abs();
        assert!(dp > 0.5, "witness should separate clearly, got Δp = {dp}");
    }

    /// Sampled XEB (the actual experimental protocol) also cannot
    /// separate the two devices beyond shot noise.
    #[test]
    fn sampled_xeb_matches_within_shot_noise() {
        let c = random_qv_circuit(5, 5, 3);
        let ideal = NativeStateVectorBackend;
        let wrong = ConjugatedStateVectorBackend::new(ConjugationScope::AllGates);
        let f_correct = xeb_score_sampled(&c, &ideal, &ideal, SHOTS, 99).unwrap();
        let f_wrong = xeb_score_sampled(&c, &ideal, &wrong, SHOTS, 99).unwrap();
        // Same seed + identical distribution ⇒ identical samples ⇒
        // identical score, exactly.
        assert!((f_correct - f_wrong).abs() < 1e-12);
    }

    /// HOG (Quantum Volume) is blind too.
    #[test]
    fn heavy_output_probability_is_blind_to_conjugation() {
        let c = random_qv_circuit(5, 5, 11);
        let ideal = NativeStateVectorBackend;
        let wrong = ConjugatedStateVectorBackend::new(ConjugationScope::AllGates);
        let h_correct = heavy_output_probability(&c, &ideal, &ideal, SHOTS, 5).unwrap();
        let h_wrong = heavy_output_probability(&c, &ideal, &wrong, SHOTS, 5).unwrap();
        assert!((h_correct - h_wrong).abs() < 1e-12);
        assert!(h_correct > 0.6, "HOG of ideal device should be heavy-biased");
    }
}
