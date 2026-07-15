//! Mirror-circuit benchmarking (Proctor et al., Sandia) and exact
//! circuit inversion over the CleitonForge IR.
//!
//! A mirror benchmark runs C followed by C† and measures the survival
//! probability P(|0…0⟩) = |⟨0|C†C|0⟩|². Since (C*)†C* = (C†C)*, a
//! simulator with a consistent conjugation bug returns exactly the same
//! survival probability as a correct one — mirror circuits are blind to
//! conjugation, like QV and XEB.

use cforge_core::{Circuit, GateKind, Operation};

use crate::trait_def::{BackendError, SimulationBackend};

/// Returns the exact inverse C† of a circuit: operations reversed, each
/// gate replaced by its adjoint.
///
/// Adjoints stay within the OpenQASM 3 `stdgates.inc` set, with one
/// exception: Csx† (controlled-√X-dagger) has no dedicated kind and is
/// emitted as Csx repeated three times (Sx⁴ = I, so Sx³ = Sx†).
pub fn inverse_circuit(circuit: &Circuit) -> Circuit {
    let mut inv = Circuit::new(circuit.num_qubits());
    for op in circuit.operations.iter().rev() {
        match op.kind {
            // Self-inverse gates.
            GateKind::Id
            | GateKind::X
            | GateKind::Y
            | GateKind::Z
            | GateKind::H
            | GateKind::Cx
            | GateKind::Cy
            | GateKind::Cz
            | GateKind::Ch
            | GateKind::Swap
            | GateKind::Ccx
            | GateKind::Cswap => {
                inv.push(Operation::new(op.kind, op.qubits.clone(), vec![]));
            }

            // Dagger pairs.
            GateKind::S => inv.push(Operation::new(GateKind::Sdg, op.qubits.clone(), vec![])),
            GateKind::Sdg => inv.push(Operation::new(GateKind::S, op.qubits.clone(), vec![])),
            GateKind::T => inv.push(Operation::new(GateKind::Tdg, op.qubits.clone(), vec![])),
            GateKind::Tdg => inv.push(Operation::new(GateKind::T, op.qubits.clone(), vec![])),
            GateKind::Sx => inv.push(Operation::new(GateKind::Sxdg, op.qubits.clone(), vec![])),
            GateKind::Sxdg => inv.push(Operation::new(GateKind::Sx, op.qubits.clone(), vec![])),

            // Rotations invert by negating the angle.
            GateKind::Rx
            | GateKind::Ry
            | GateKind::Rz
            | GateKind::Phase
            | GateKind::Crx
            | GateKind::Cry
            | GateKind::Crz
            | GateKind::Cp => {
                inv.push(Operation::new(
                    op.kind,
                    op.qubits.clone(),
                    vec![-op.params[0]],
                ));
            }

            // U(θ,φ,λ)† = U(−θ,−λ,−φ).
            GateKind::U => inv.push(Operation::new(
                GateKind::U,
                op.qubits.clone(),
                vec![-op.params[0], -op.params[2], -op.params[1]],
            )),
            // CU(θ,φ,λ,γ)† = CU(−θ,−λ,−φ,−γ).
            GateKind::Cu => inv.push(Operation::new(
                GateKind::Cu,
                op.qubits.clone(),
                vec![-op.params[0], -op.params[2], -op.params[1], -op.params[3]],
            )),

            // Csx† = Csx³.
            GateKind::Csx => {
                for _ in 0..3 {
                    inv.push(Operation::new(GateKind::Csx, op.qubits.clone(), vec![]));
                }
            }
        }
    }
    inv
}

/// Builds the mirror circuit C†C and returns the survival probability
/// P(|0…0⟩) as computed by `backend`. A correct backend returns 1.0
/// exactly (up to f64 accumulation); so does a consistently conjugated
/// one, which is precisely the blindness.
pub fn mirror_survival(
    circuit: &Circuit,
    backend: &dyn SimulationBackend,
) -> Result<f64, BackendError> {
    let mut mirrored = circuit.clone();
    mirrored
        .operations
        .extend(inverse_circuit(circuit).operations);
    let result = backend.run(&mirrored, 0, 0)?;
    Ok(result.statevector[0].norm_sqr())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::statevector::NativeStateVectorBackend;
    use crate::xeb::{random_qv_circuit, ConjugatedStateVectorBackend, ConjugationScope};

    /// C†C = I for a circuit exercising every gate kind.
    #[test]
    fn inverse_undoes_every_gate_kind() {
        let mut c = Circuit::new(3);
        let ops: Vec<(GateKind, Vec<usize>, Vec<f64>)> = vec![
            (GateKind::H, vec![0], vec![]),
            (GateKind::X, vec![1], vec![]),
            (GateKind::Y, vec![2], vec![]),
            (GateKind::Z, vec![0], vec![]),
            (GateKind::S, vec![1], vec![]),
            (GateKind::Sdg, vec![2], vec![]),
            (GateKind::T, vec![0], vec![]),
            (GateKind::Tdg, vec![1], vec![]),
            (GateKind::Sx, vec![2], vec![]),
            (GateKind::Sxdg, vec![0], vec![]),
            (GateKind::Rx, vec![1], vec![0.3]),
            (GateKind::Ry, vec![2], vec![0.7]),
            (GateKind::Rz, vec![0], vec![1.1]),
            (GateKind::Phase, vec![1], vec![0.5]),
            (GateKind::U, vec![2], vec![0.4, 0.8, 1.2]),
            (GateKind::Cx, vec![0, 1], vec![]),
            (GateKind::Cy, vec![1, 2], vec![]),
            (GateKind::Cz, vec![2, 0], vec![]),
            (GateKind::Ch, vec![0, 2], vec![]),
            (GateKind::Csx, vec![1, 0], vec![]),
            (GateKind::Crx, vec![2, 1], vec![0.6]),
            (GateKind::Cry, vec![0, 1], vec![0.9]),
            (GateKind::Crz, vec![1, 2], vec![1.3]),
            (GateKind::Cp, vec![2, 0], vec![0.2]),
            (GateKind::Cu, vec![0, 1], vec![0.5, 0.6, 0.7, 0.8]),
            (GateKind::Swap, vec![1, 2], vec![]),
            (GateKind::Ccx, vec![0, 1, 2], vec![]),
            (GateKind::Cswap, vec![2, 0, 1], vec![]),
        ];
        for (kind, qubits, params) in ops {
            c.push(Operation::new(kind, qubits, params));
        }
        let survival = mirror_survival(&c, &NativeStateVectorBackend).unwrap();
        assert!(
            (survival - 1.0).abs() < 1e-9,
            "C†C should return to |000⟩, got survival {survival}"
        );
    }

    /// Corollary (mirror blindness): the conjugated simulator also
    /// reports perfect survival — the benchmark cannot see the bug.
    #[test]
    fn mirror_is_blind_to_conjugation() {
        let c = random_qv_circuit(4, 4, 21);
        let good = mirror_survival(&c, &NativeStateVectorBackend).unwrap();
        let bad = mirror_survival(
            &c,
            &ConjugatedStateVectorBackend::new(ConjugationScope::AllGates),
        )
        .unwrap();
        assert!((good - 1.0).abs() < 1e-9);
        assert!((bad - 1.0).abs() < 1e-9);
    }
}
