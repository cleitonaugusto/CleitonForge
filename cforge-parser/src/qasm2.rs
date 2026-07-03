//! OpenQASM 2.0 parser wrapping the `qasm` crate.
//!
//! Translates a QASM 2 source string into CleitonForge's canonical IR. Only
//! gate applications and qubit register declarations are translated;
//! classical registers, measurements, barriers, and conditionals are
//! silently skipped (they have no bearing on the unitary part of the
//! circuit that backends simulate).

use std::collections::HashMap;
use std::path::Path;

use cforge_core::{Circuit, Operation};

use crate::error::ParseError;
use crate::translate::translate_gate;

/// Parses an OpenQASM 2.0 source string and returns the canonical circuit.
///
/// `search_dir` is the directory used to resolve `include` statements (pass
/// `std::env::current_dir().unwrap()` when there are none).
pub fn parse_qasm2(source: &str, search_dir: &Path) -> Result<Circuit, ParseError> {
    let processed = qasm::process(source, search_dir);
    let mut tokens = qasm::lex(&processed);
    let ast = qasm::parse(&mut tokens)
        .map_err(|e| ParseError::SyntaxError(format!("{e:?}")))?;

    // First pass: collect qubit register declarations and assign flat indices.
    // reg_map: register_name → (base_flat_index, size)
    let mut reg_map: HashMap<String, (usize, usize)> = HashMap::new();
    let mut total_qubits: usize = 0;

    for node in &ast {
        if let qasm::AstNode::QReg(name, size) = node {
            let size = *size as usize;
            reg_map.insert(name.clone(), (total_qubits, size));
            total_qubits += size;
        }
    }

    let mut circuit = Circuit::new(total_qubits);

    // Second pass: translate gate applications.
    for node in &ast {
        let qasm::AstNode::ApplyGate(gate_name, args, params) = node else {
            continue;
        };

        // Evaluate params first so translate_gate can rewrite them (e.g. u2).
        let raw_params = params
            .iter()
            .map(|s| eval_param(s))
            .collect::<Result<Vec<_>, _>>()?;

        let (kind, params) = translate_gate(gate_name, raw_params)
            .ok_or_else(|| ParseError::UnknownGate(gate_name.clone()))?;

        // Each argument is either a single indexed qubit or a whole register.
        // Whole-register args broadcast: `h q;` on q[3] → h q[0]; h q[1]; h q[2].
        let arg_lists: Vec<Vec<usize>> = args
            .iter()
            .map(|arg| resolve_arg(arg, &reg_map))
            .collect::<Result<_, _>>()?;

        // broadcast_n: 1 for all-indexed, N for register application.
        let broadcast_n = arg_lists.iter().map(|v| v.len()).max().unwrap_or(1);

        // All register-typed args must agree on size.
        for v in &arg_lists {
            if v.len() != 1 && v.len() != broadcast_n {
                return Err(ParseError::SyntaxError(format!(
                    "register size mismatch in gate '{gate_name}'"
                )));
            }
        }

        for i in 0..broadcast_n {
            let qubits: Vec<usize> = arg_lists
                .iter()
                .map(|v| if v.len() == 1 { v[0] } else { v[i] })
                .collect();
            let mut op = Operation::new(kind, qubits, params.clone());
            op.source_framework = Some("openqasm2".to_string());
            op.original_gate_name = Some(gate_name.clone());
            circuit.push(op);
        }
    }

    circuit.validate().map_err(|e| ParseError::SyntaxError(format!("{e:?}")))?;
    Ok(circuit)
}

/// Resolves a gate argument to a list of flat qubit indices.
///
/// - `Qubit(name, i)` → `[base + i]` (single qubit)
/// - `Register(name)` → `[base, base+1, …, base+size-1]` (whole-register
///   broadcast: `h q;` on a 3-qubit register expands to three H gates)
fn resolve_arg(
    arg: &qasm::Argument,
    reg_map: &HashMap<String, (usize, usize)>,
) -> Result<Vec<usize>, ParseError> {
    match arg {
        qasm::Argument::Qubit(name, index) => {
            let (base, size) = reg_map
                .get(name)
                .ok_or_else(|| ParseError::UndeclaredQubit(name.clone()))?;
            let idx = *index as usize;
            if idx >= *size {
                return Err(ParseError::UndeclaredQubit(format!("{name}[{idx}]")));
            }
            Ok(vec![base + idx])
        }
        qasm::Argument::Register(name) => {
            let (base, size) = reg_map
                .get(name)
                .ok_or_else(|| ParseError::UndeclaredQubit(name.clone()))?;
            Ok((*base..*base + *size).collect())
        }
    }
}

/// Evaluates a QASM 2 parameter expression string to an `f64`.
///
/// Handles: numeric literals, `pi`, and the binary operators +/-/*//
/// with `pi`. Complex expressions beyond these patterns are rejected.
fn eval_param(expr: &str) -> Result<f64, ParseError> {
    let s = expr.trim();
    // Try plain float first.
    if let Ok(v) = s.parse::<f64>() {
        return Ok(v);
    }
    // Simple substitution-then-eval for expressions involving `pi`.
    let substituted = s.replace("pi", &std::f64::consts::PI.to_string());
    if let Ok(v) = eval_simple_expr(&substituted) {
        return Ok(v);
    }
    Err(ParseError::InvalidParam(expr.to_string()))
}

/// Evaluates simple additive/multiplicative expressions over f64 literals.
fn eval_simple_expr(s: &str) -> Result<f64, ()> {
    // Handle +/-
    for (i, c) in s.char_indices().rev() {
        if (c == '+' || c == '-') && i > 0 {
            let left = eval_simple_expr(s[..i].trim())?;
            let right: f64 = eval_simple_expr(s[i + 1..].trim())?;
            return Ok(if c == '+' { left + right } else { left - right });
        }
    }
    // Handle * and /
    for (i, c) in s.char_indices().rev() {
        if c == '*' || c == '/' {
            let left = eval_simple_expr(s[..i].trim())?;
            let right = eval_simple_expr(s[i + 1..].trim())?;
            return Ok(if c == '*' { left * right } else { left / right });
        }
    }
    s.parse::<f64>().map_err(|_| ())
}

#[cfg(test)]
mod tests {
    use super::*;
    use cforge_core::GateKind;

    // No include — avoids requiring qelib1.inc on disk during tests.
    const BELL_QASM2: &str = r#"
OPENQASM 2.0;
qreg q[2];
creg c[2];
h q[0];
cx q[0],q[1];
measure q[0] -> c[0];
measure q[1] -> c[1];
"#;

    #[test]
    fn bell_state_qasm2() {
        let dir = std::env::current_dir().unwrap();
        let circuit = parse_qasm2(BELL_QASM2, &dir).unwrap();
        assert_eq!(circuit.num_qubits(), 2);
        assert_eq!(circuit.operations.len(), 2);
        assert_eq!(circuit.operations[0].kind, GateKind::H);
        assert_eq!(circuit.operations[0].qubits, vec![0]);
        assert_eq!(circuit.operations[1].kind, GateKind::Cx);
        assert_eq!(circuit.operations[1].qubits, vec![0, 1]);
    }

    #[test]
    fn parameterized_gate_qasm2() {
        let source = "OPENQASM 2.0;\nqreg q[1];\nrz(pi/2) q[0];\n";
        let dir = std::env::current_dir().unwrap();
        let circuit = parse_qasm2(source, &dir).unwrap();
        assert_eq!(circuit.operations.len(), 1);
        assert_eq!(circuit.operations[0].kind, GateKind::Rz);
        let expected = std::f64::consts::FRAC_PI_2;
        assert!((circuit.operations[0].params[0] - expected).abs() < 1e-10);
    }

    #[test]
    fn whole_register_single_qubit_gate() {
        // `h q;` on a 3-qubit register must expand to h q[0]; h q[1]; h q[2].
        let source = "OPENQASM 2.0;\nqreg q[3];\nh q;\n";
        let dir = std::env::current_dir().unwrap();
        let circuit = parse_qasm2(source, &dir).unwrap();
        assert_eq!(circuit.num_qubits(), 3);
        assert_eq!(circuit.operations.len(), 3);
        for (i, op) in circuit.operations.iter().enumerate() {
            assert_eq!(op.kind, GateKind::H);
            assert_eq!(op.qubits, vec![i]);
        }
    }

    #[test]
    fn whole_register_two_qubit_gate() {
        // `cx q, r;` where both are 2-qubit registers →
        // cx q[0],r[0]; cx q[1],r[1];
        let source = "OPENQASM 2.0;\nqreg q[2];\nqreg r[2];\ncx q,r;\n";
        let dir = std::env::current_dir().unwrap();
        let circuit = parse_qasm2(source, &dir).unwrap();
        assert_eq!(circuit.num_qubits(), 4);
        assert_eq!(circuit.operations.len(), 2);
        assert_eq!(circuit.operations[0].qubits, vec![0, 2]); // q[0], r[0]
        assert_eq!(circuit.operations[1].qubits, vec![1, 3]); // q[1], r[1]
    }

    #[test]
    fn u2_gate_rewrites_to_u_qasm2() {
        // u2(phi, lam) is common in IBM/Qiskit exported circuits.
        // It must parse as U(pi/2, phi, lam) rather than UnknownGate.
        let source = "OPENQASM 2.0;\nqreg q[1];\nu2(0,pi) q[0];\n";
        let dir = std::env::current_dir().unwrap();
        let circuit = parse_qasm2(source, &dir).unwrap();
        assert_eq!(circuit.operations.len(), 1);
        let op = &circuit.operations[0];
        assert_eq!(op.kind, GateKind::U);
        assert_eq!(op.params.len(), 3);
        assert!((op.params[0] - std::f64::consts::FRAC_PI_2).abs() < 1e-10);
        assert!((op.params[1] - 0.0).abs() < 1e-10);
        assert!((op.params[2] - std::f64::consts::PI).abs() < 1e-10);
    }
}
