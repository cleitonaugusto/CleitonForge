//! Typed random circuit generation over the cforge IR.
//!
//! Convention-sensitive gates (the diagonal phase family and the
//! discrete complex gates) are weighted 5× relative to plain Paulis:
//! phase-convention bugs only surface when a bugged gate meets a
//! correctly-implemented complex gate, so the generator biases the
//! search toward exactly those mixtures.

use std::f64::consts::PI;

use cforge_core::{Circuit, GateKind, Operation};
use rand::rngs::StdRng;
use rand::Rng;

/// Gate pool with sampling weights. Phase-family and complex discrete
/// gates carry weight 5; structural/real gates fill out coverage.
const GATE_WEIGHTS: &[(GateKind, u32)] = &[
    // Phase family — the primary bug habitat.
    (GateKind::Rz, 5),
    (GateKind::Phase, 5),
    (GateKind::Crz, 5),
    (GateKind::Cp, 5),
    // Discrete complex gates — the "phase reference" that makes partial
    // bugs visible in probabilities.
    (GateKind::S, 5),
    (GateKind::Sdg, 5),
    (GateKind::T, 5),
    (GateKind::Tdg, 5),
    (GateKind::Sx, 5),
    (GateKind::Sxdg, 5),
    // Superposition makers and rotations.
    (GateKind::H, 4),
    (GateKind::Rx, 2),
    (GateKind::Ry, 2),
    (GateKind::U, 2),
    // Entanglers.
    (GateKind::Cx, 3),
    (GateKind::Cy, 1),
    (GateKind::Cz, 1),
    (GateKind::Ch, 1),
    (GateKind::Csx, 2),
    (GateKind::Crx, 1),
    (GateKind::Cry, 1),
    (GateKind::Cu, 1),
    // Real structural gates.
    (GateKind::X, 1),
    (GateKind::Y, 1),
    (GateKind::Z, 1),
    (GateKind::Swap, 1),
    (GateKind::Ccx, 1),
    (GateKind::Cswap, 1),
];

fn pick_gate(rng: &mut StdRng, max_qubits: usize) -> GateKind {
    let eligible: Vec<&(GateKind, u32)> = GATE_WEIGHTS
        .iter()
        .filter(|(kind, _)| kind.num_qubits() <= max_qubits)
        .collect();
    let total: u32 = eligible.iter().map(|(_, w)| w).sum();
    let mut roll = rng.gen_range(0..total);
    for (kind, weight) in &eligible {
        if roll < *weight {
            return *kind;
        }
        roll -= weight;
    }
    eligible.last().unwrap().0
}

/// Generates a random circuit of exactly `depth` gates over
/// `num_qubits` qubits, deterministic in the RNG state.
pub fn random_circuit(rng: &mut StdRng, num_qubits: usize, depth: usize) -> Circuit {
    assert!(num_qubits > 0, "need at least one qubit");
    let mut c = Circuit::new(num_qubits);
    for _ in 0..depth {
        let gate = pick_gate(rng, num_qubits);
        let nq = gate.num_qubits();

        // Distinct qubit operands.
        let mut available: Vec<usize> = (0..num_qubits).collect();
        let mut qubits = Vec::with_capacity(nq);
        for _ in 0..nq {
            let idx = rng.gen_range(0..available.len());
            qubits.push(available.remove(idx));
        }

        let params: Vec<f64> = (0..gate.num_params())
            .map(|_| rng.gen_range(-PI..PI))
            .collect();

        c.push(Operation::new(gate, qubits, params));
    }
    c
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    #[test]
    fn generated_circuits_are_valid() {
        let mut rng = StdRng::seed_from_u64(5);
        for _ in 0..500 {
            let c = random_circuit(&mut rng, 3, 10);
            assert_eq!(c.validate(), Ok(()));
            assert_eq!(c.operations.len(), 10);
        }
    }

    #[test]
    fn single_qubit_circuits_only_use_single_qubit_gates() {
        let mut rng = StdRng::seed_from_u64(9);
        for _ in 0..200 {
            let c = random_circuit(&mut rng, 1, 6);
            assert_eq!(c.validate(), Ok(()));
        }
    }

    #[test]
    fn generation_is_deterministic_in_seed() {
        let a = random_circuit(&mut StdRng::seed_from_u64(77), 3, 8);
        let b = random_circuit(&mut StdRng::seed_from_u64(77), 3, 8);
        assert_eq!(a.operations, b.operations);
    }
}
