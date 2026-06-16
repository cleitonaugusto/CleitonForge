//! Canonical circuit representation shared by every parser, backend, and
//! metric in CleitonForge. Parsers translate framework-specific circuits
//! into this IR; backends and metrics only ever see this IR, never a
//! framework-native representation, which is what makes results
//! comparable across frameworks.

use crate::gate::GateKind;

/// A single qubit in a [`Circuit`], identified by its index.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Qubit {
    pub index: usize,
    /// Optional name as it appeared in the source (e.g. `q[0]`).
    pub name: Option<String>,
}

/// A single gate application within a [`Circuit`].
#[derive(Debug, Clone, PartialEq)]
pub struct Operation {
    pub kind: GateKind,
    /// Indices of the qubits this operation acts on, in the order
    /// expected by `kind` (e.g. for `Cx`, `qubits[0]` is the control and
    /// `qubits[1]` is the target).
    pub qubits: Vec<usize>,
    /// Numeric parameters, e.g. rotation angles.
    pub params: Vec<f64>,
    /// Name of the framework this operation was translated from, if known.
    pub source_framework: Option<String>,
    /// The gate name as it appeared in the source framework, if it
    /// differs from the canonical `qasm_name`.
    pub original_gate_name: Option<String>,
}

impl Operation {
    /// Builds an operation with no provenance metadata.
    pub fn new(kind: GateKind, qubits: Vec<usize>, params: Vec<f64>) -> Self {
        Self {
            kind,
            qubits,
            params,
            source_framework: None,
            original_gate_name: None,
        }
    }
}

/// An error describing why a [`Circuit`] is not well-formed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CircuitError {
    /// An operation referenced a qubit index outside `0..num_qubits`.
    QubitOutOfRange { operation_index: usize, qubit: usize },
    /// An operation's qubit count didn't match what its `GateKind` expects.
    WrongQubitCount {
        operation_index: usize,
        expected: usize,
        found: usize,
    },
    /// An operation's parameter count didn't match what its `GateKind` expects.
    WrongParamCount {
        operation_index: usize,
        expected: usize,
        found: usize,
    },
}

/// A quantum circuit in CleitonForge's canonical intermediate
/// representation: a fixed set of qubits and a sequence of gate
/// applications over them.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Circuit {
    pub qubits: Vec<Qubit>,
    pub operations: Vec<Operation>,
}

impl Circuit {
    /// Creates an empty circuit with `num_qubits` unnamed qubits.
    pub fn new(num_qubits: usize) -> Self {
        Self {
            qubits: (0..num_qubits).map(|index| Qubit { index, name: None }).collect(),
            operations: Vec::new(),
        }
    }

    pub fn num_qubits(&self) -> usize {
        self.qubits.len()
    }

    /// Appends an operation to the circuit.
    pub fn push(&mut self, operation: Operation) {
        self.operations.push(operation);
    }

    /// Checks that every operation references valid qubit indices and has
    /// the qubit/parameter counts its `GateKind` expects.
    pub fn validate(&self) -> Result<(), CircuitError> {
        let num_qubits = self.num_qubits();
        for (operation_index, operation) in self.operations.iter().enumerate() {
            for &qubit in &operation.qubits {
                if qubit >= num_qubits {
                    return Err(CircuitError::QubitOutOfRange { operation_index, qubit });
                }
            }

            let expected_qubits = operation.kind.num_qubits();
            if operation.qubits.len() != expected_qubits {
                return Err(CircuitError::WrongQubitCount {
                    operation_index,
                    expected: expected_qubits,
                    found: operation.qubits.len(),
                });
            }

            let expected_params = operation.kind.num_params();
            if operation.params.len() != expected_params {
                return Err(CircuitError::WrongParamCount {
                    operation_index,
                    expected: expected_params,
                    found: operation.params.len(),
                });
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bell_state_circuit() -> Circuit {
        let mut circuit = Circuit::new(2);
        circuit.push(Operation::new(GateKind::H, vec![0], vec![]));
        circuit.push(Operation::new(GateKind::Cx, vec![0, 1], vec![]));
        circuit
    }

    #[test]
    fn bell_state_has_expected_shape() {
        let circuit = bell_state_circuit();
        assert_eq!(circuit.num_qubits(), 2);
        assert_eq!(circuit.operations.len(), 2);
        assert_eq!(circuit.operations[0].kind, GateKind::H);
        assert_eq!(circuit.operations[1].kind, GateKind::Cx);
    }

    #[test]
    fn bell_state_is_valid() {
        assert_eq!(bell_state_circuit().validate(), Ok(()));
    }

    #[test]
    fn rejects_out_of_range_qubit() {
        let mut circuit = Circuit::new(1);
        circuit.push(Operation::new(GateKind::Cx, vec![0, 1], vec![]));
        assert_eq!(
            circuit.validate(),
            Err(CircuitError::QubitOutOfRange { operation_index: 0, qubit: 1 })
        );
    }

    #[test]
    fn rejects_wrong_qubit_count() {
        let mut circuit = Circuit::new(2);
        circuit.push(Operation::new(GateKind::H, vec![0, 1], vec![]));
        assert_eq!(
            circuit.validate(),
            Err(CircuitError::WrongQubitCount { operation_index: 0, expected: 1, found: 2 })
        );
    }

    #[test]
    fn rejects_wrong_param_count() {
        let mut circuit = Circuit::new(1);
        circuit.push(Operation::new(GateKind::Rz, vec![0], vec![]));
        assert_eq!(
            circuit.validate(),
            Err(CircuitError::WrongParamCount { operation_index: 0, expected: 1, found: 0 })
        );
    }

    #[test]
    fn provenance_metadata_is_optional_but_settable() {
        let mut operation = Operation::new(GateKind::X, vec![0], vec![]);
        operation.source_framework = Some("qiskit".to_string());
        operation.original_gate_name = Some("x".to_string());
        assert_eq!(operation.source_framework.as_deref(), Some("qiskit"));
    }
}
