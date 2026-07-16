//! cforge-fuzz — phase-sensitive differential fuzzing for quantum
//! simulators.
//!
//! The library implements the four stages of a fuzzing campaign:
//!
//! 1. [`generator`] — typed random circuits over the cforge IR, with
//!    convention-sensitive gates weighted 5×;
//! 2. [`oracle`] — a strict three-level hierarchy: amplitude equality
//!    modulo global phase (N1) ⊃ probability equality (N2) ⊃ Z-observable
//!    equality (N3). Bugs at N1 but not N2 are *provably invisible* to
//!    every sampling benchmark (QV, XEB, mirror, RB, QPE) by the
//!    conjugation-invariance theorem;
//! 3. [`shrinker`] — delta-debugging to a minimal counterexample;
//! 4. [`triage`] — classification by QGCS dimension and benchmark
//!    visibility, plus [`zoo`] for the public bug corpus.

pub mod generator;
pub mod oracle;
pub mod shrinker;
pub mod triage;
pub mod zoo;

#[cfg(test)]
mod tests {
    use crate::generator::random_circuit;
    use crate::oracle::{distance, OracleLevel};
    use crate::shrinker::shrink;
    use cforge_backends::{
        ConjugatedStateVectorBackend, ConjugationScope, NativeStateVectorBackend,
        SimulationBackend,
    };
    use cforge_core::{Circuit, GateKind};
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    const TOL: f64 = 1e-6;

    fn diverges(c: &Circuit, device: &dyn SimulationBackend, level: OracleLevel) -> Option<f64> {
        let a = NativeStateVectorBackend.run(c, 0, 0).ok()?.statevector;
        let b = device.run(c, 0, 0).ok()?.statevector;
        let d = distance(&a, &b, level);
        (d > TOL).then_some(d)
    }

    /// The theorem, live: with the probability oracle (what sampling
    /// benchmarks measure), the fully conjugated simulator survives
    /// thousands of random circuits without a single divergence.
    #[test]
    fn probability_oracle_is_blind_to_global_conjugation() {
        let device = ConjugatedStateVectorBackend::new(ConjugationScope::AllGates);
        let mut rng = StdRng::seed_from_u64(1);
        for _ in 0..2_000 {
            let c = random_circuit(&mut rng, 3, 8);
            assert!(
                diverges(&c, &device, OracleLevel::Probability).is_none(),
                "theorem violated?! circuit: {:?}",
                c.operations
            );
        }
    }

    /// The amplitude oracle separates the same device immediately, and
    /// the shrinker reduces the witness to ≤ 3 gates.
    #[test]
    fn amplitude_oracle_finds_global_conjugation_fast() {
        let device = ConjugatedStateVectorBackend::new(ConjugationScope::AllGates);
        let mut rng = StdRng::seed_from_u64(1);
        let mut tried = 0;
        for _ in 0..500 {
            tried += 1;
            let c = random_circuit(&mut rng, 3, 8);
            if diverges(&c, &device, OracleLevel::Amplitude).is_some() {
                let minimal = shrink(&c, &|cand| {
                    diverges(cand, &device, OracleLevel::Amplitude).is_some()
                });
                assert!(
                    minimal.operations.len() <= 3,
                    "witness should shrink to ≤3 gates, got {}",
                    minimal.operations.len()
                );
                return;
            }
        }
        panic!("no divergence found in {tried} circuits");
    }

    /// The partial (quantrs2-class) bug IS visible to the probability
    /// oracle — but only via circuits mixing the bugged rotation with a
    /// correct complex gate. The fuzzer must rediscover that witness
    /// structure on its own.
    #[test]
    fn probability_oracle_finds_partial_bug_and_witness_mixes_gates() {
        let device = ConjugatedStateVectorBackend::new(ConjugationScope::RotationSignOnly);
        let mut rng = StdRng::seed_from_u64(7);
        for _ in 0..20_000 {
            let c = random_circuit(&mut rng, 3, 8);
            if diverges(&c, &device, OracleLevel::Probability).is_some() {
                let minimal = shrink(&c, &|cand| {
                    diverges(cand, &device, OracleLevel::Probability).is_some()
                });
                let has_bugged_rotation = minimal.operations.iter().any(|op| {
                    matches!(
                        op.kind,
                        GateKind::Rz | GateKind::Crz | GateKind::Phase | GateKind::Cp
                    )
                });
                assert!(
                    has_bugged_rotation,
                    "minimal witness must exercise the bugged family: {:?}",
                    minimal.operations
                );
                assert!(
                    minimal.operations.len() <= 6,
                    "witness should be small, got {} gates",
                    minimal.operations.len()
                );
                return;
            }
        }
        panic!("partial bug not found — generator or oracle regressed");
    }

    /// Global phase must NOT count as divergence at any level: a correct
    /// device that differs only by e^{iφ} is not a bug.
    #[test]
    fn global_phase_is_not_a_divergence() {
        use num_complex::Complex64;
        let a = vec![
            Complex64::new(0.6, 0.0),
            Complex64::new(0.0, 0.8),
        ];
        let phase = Complex64::from_polar(1.0, 1.234);
        let b: Vec<Complex64> = a.iter().map(|x| x * phase).collect();
        assert!(distance(&a, &b, OracleLevel::Amplitude) < 1e-12);
        assert!(distance(&a, &b, OracleLevel::Probability) < 1e-12);
    }
}
