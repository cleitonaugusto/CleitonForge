//! Delta-debugging shrinker: reduces a diverging circuit to a minimal
//! counterexample while preserving the failure.
//!
//! Two passes: chunked removal (ddmin-style halving) for speed on long
//! circuits, then greedy single-gate removal until a fixed point. The
//! result is 1-minimal — removing any single remaining gate makes the
//! divergence disappear — which is what makes a bug report convincing.

use cforge_core::{Circuit, Operation};

fn rebuild(num_qubits: usize, ops: &[Operation]) -> Circuit {
    let mut c = Circuit::new(num_qubits);
    for op in ops {
        c.push(op.clone());
    }
    c
}

/// Shrinks `circuit` to a 1-minimal counterexample under `still_fails`.
///
/// `still_fails` must return true for the input circuit; the returned
/// circuit also satisfies it.
pub fn shrink(circuit: &Circuit, still_fails: &dyn Fn(&Circuit) -> bool) -> Circuit {
    let nq = circuit.num_qubits();
    let mut current: Vec<Operation> = circuit.operations.clone();

    // Pass 1: remove contiguous chunks, halving the chunk size.
    let mut chunk = current.len() / 2;
    while chunk >= 1 {
        let mut start = 0;
        while start < current.len() {
            let end = (start + chunk).min(current.len());
            let mut candidate = current.clone();
            candidate.drain(start..end);
            if !candidate.is_empty() && still_fails(&rebuild(nq, &candidate)) {
                current = candidate;
                // Re-scan from the beginning at this chunk size.
                start = 0;
            } else {
                start += chunk;
            }
        }
        chunk /= 2;
    }

    // Pass 2: greedy single-gate removal to a fixed point (1-minimal).
    let mut changed = true;
    while changed {
        changed = false;
        let mut i = 0;
        while i < current.len() {
            let mut candidate = current.clone();
            candidate.remove(i);
            if !candidate.is_empty() && still_fails(&rebuild(nq, &candidate)) {
                current = candidate;
                changed = true;
            } else {
                i += 1;
            }
        }
    }

    rebuild(nq, &current)
}

#[cfg(test)]
mod tests {
    use super::*;
    use cforge_core::GateKind;

    /// Failure predicate: "circuit contains both an Rz and an S".
    /// The minimal circuit under it has exactly 2 gates.
    #[test]
    fn shrinks_to_the_minimal_failing_core() {
        let mut c = Circuit::new(2);
        for kind in [
            GateKind::H,
            GateKind::X,
            GateKind::Rz,
            GateKind::Y,
            GateKind::Z,
            GateKind::S,
            GateKind::H,
            GateKind::X,
        ] {
            let params = if kind == GateKind::Rz { vec![0.5] } else { vec![] };
            c.push(Operation::new(kind, vec![0], params));
        }

        let fails = |cand: &Circuit| {
            let has_rz = cand.operations.iter().any(|o| o.kind == GateKind::Rz);
            let has_s = cand.operations.iter().any(|o| o.kind == GateKind::S);
            has_rz && has_s
        };
        let minimal = shrink(&c, &fails);
        assert_eq!(minimal.operations.len(), 2);
        assert!(fails(&minimal));
    }

    /// The shrinker must return a circuit that still fails, never an
    /// empty or passing one.
    #[test]
    fn result_still_fails() {
        let mut c = Circuit::new(1);
        for _ in 0..7 {
            c.push(Operation::new(GateKind::H, vec![0], vec![]));
        }
        let fails = |cand: &Circuit| cand.operations.len() >= 3;
        let minimal = shrink(&c, &fails);
        assert_eq!(minimal.operations.len(), 3);
    }
}
