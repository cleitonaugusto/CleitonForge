//! Circuit depth and gate-count metrics, always computed on the
//! canonical IR so results are comparable across backends.

use std::collections::HashMap;
use cforge_core::{Circuit, GateKind};

/// Per-gate-kind count and total, plus circuit depth.
#[derive(Debug, Clone, Default)]
pub struct CircuitStats {
    /// Total number of gate operations.
    pub gate_count: usize,
    /// Circuit depth: number of sequential layers when independent
    /// operations on different qubits are parallelized.
    pub depth: usize,
    /// How many times each `GateKind` appears.
    pub gate_counts_by_kind: HashMap<GateKind, usize>,
}

/// Computes circuit statistics from the canonical IR.
///
/// Depth is calculated as the critical path: for each operation, its
/// layer is 1 + max(layer of previous operations on any of its qubits).
pub fn compute_stats(circuit: &Circuit) -> CircuitStats {
    let n = circuit.num_qubits();
    // Layer at which each qubit is last occupied.
    let mut qubit_layer = vec![0usize; n];

    let mut gate_counts_by_kind: HashMap<GateKind, usize> = HashMap::new();
    let mut max_depth = 0usize;

    for op in &circuit.operations {
        let op_layer = op.qubits.iter().map(|&q| qubit_layer[q]).max().unwrap_or(0) + 1;
        for &q in &op.qubits {
            qubit_layer[q] = op_layer;
        }
        max_depth = max_depth.max(op_layer);
        *gate_counts_by_kind.entry(op.kind).or_insert(0) += 1;
    }

    CircuitStats {
        gate_count: circuit.operations.len(),
        depth: max_depth,
        gate_counts_by_kind,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cforge_core::Operation;

    fn bell() -> Circuit {
        let mut c = Circuit::new(2);
        c.push(Operation::new(GateKind::H, vec![0], vec![]));
        c.push(Operation::new(GateKind::Cx, vec![0, 1], vec![]));
        c
    }

    #[test]
    fn bell_state_gate_count_and_depth() {
        let stats = compute_stats(&bell());
        assert_eq!(stats.gate_count, 2);
        // H on q[0] (layer 1), CX on q[0]+q[1] (layer 2, must wait for H)
        assert_eq!(stats.depth, 2);
    }

    #[test]
    fn parallel_gates_reduce_depth() {
        // H q[0] and H q[1] can run in parallel → depth 1
        let mut c = Circuit::new(2);
        c.push(Operation::new(GateKind::H, vec![0], vec![]));
        c.push(Operation::new(GateKind::H, vec![1], vec![]));
        let stats = compute_stats(&c);
        assert_eq!(stats.gate_count, 2);
        assert_eq!(stats.depth, 1);
    }

    #[test]
    fn gate_counts_by_kind_correct() {
        let stats = compute_stats(&bell());
        assert_eq!(stats.gate_counts_by_kind.get(&GateKind::H).copied(), Some(1));
        assert_eq!(stats.gate_counts_by_kind.get(&GateKind::Cx).copied(), Some(1));
    }
}
