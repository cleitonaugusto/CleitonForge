//! OpenQASM 3 parser wrapping the `oq3_semantics` crate.
//!
//! Translates a QASM 3 source string into CleitonForge's canonical IR.
//! Supports: qubit register declarations and gate applications from the
//! stdgates.inc gate set. Measurements, classical registers, and control
//! flow are silently skipped — they are irrelevant to the unitary
//! simulation layer.

use std::collections::HashMap;

use cforge_core::{Circuit, Operation, Qubit};
use oq3_semantics::asg::{self, Expr, Literal, BinaryOp, ArithOp, UnaryOp};
use oq3_semantics::symbols::{SymbolId, SymbolTable, SymbolType};
use oq3_semantics::syntax_to_semantics::parse_source_string;
use oq3_semantics::types::Type;

use crate::error::ParseError;
use crate::translate::gate_kind_from_name;

/// Parses an OpenQASM 3 source string and returns the canonical circuit.
pub fn parse_qasm3(source: &str) -> Result<Circuit, ParseError> {
    let result = parse_source_string(source, None, None::<&[&std::path::Path]>);

    if result.any_syntax_errors() {
        return Err(ParseError::SyntaxError(
            "OpenQASM 3 syntax error (run with verbose diagnostics for details)".to_string(),
        ));
    }

    let (program, sem_errors, symbol_table) = result.take_context().as_tuple();

    if !sem_errors.is_empty() {
        return Err(ParseError::SyntaxError(format!(
            "{} semantic error(s) in OpenQASM 3 source",
            sem_errors.len()
        )));
    }

    // First pass: collect qubit register declarations and assign flat indices.
    // qubit_base: SymbolId of the register → base flat index
    let mut qubit_base: HashMap<SymbolId, usize> = HashMap::new();
    let mut total_qubits: usize = 0;
    let mut named_qubits: Vec<Qubit> = Vec::new();

    for stmt in program.stmts() {
        let asg::Stmt::DeclareQuantum(dq) = stmt else {
            continue;
        };
        let sym_id = match dq.name() {
            Ok(id) => id.clone(),
            Err(_) => continue,
        };
        let sym = &symbol_table[&sym_id];
        let width = match sym.symbol_type() {
            Type::Qubit => 1,
            Type::QubitArray(dims) => dims.dims()[0],
            _ => continue,
        };
        qubit_base.insert(sym_id, total_qubits);
        for i in 0..width {
            named_qubits.push(Qubit {
                index: total_qubits + i,
                name: Some(format!("{}[{}]", sym.name(), i)),
            });
        }
        total_qubits += width;
    }

    let mut circuit = Circuit {
        qubits: named_qubits,
        operations: Vec::new(),
    };

    // Second pass: translate gate applications.
    for stmt in program.stmts() {
        let asg::Stmt::GateCall(gc) = stmt else {
            continue;
        };
        let sym_id = match gc.name() {
            Ok(id) => id.clone(),
            Err(_) => continue,
        };
        let gate_name = symbol_table[&sym_id].name().to_string();

        let kind = gate_kind_from_name(&gate_name)
            .ok_or_else(|| ParseError::UnknownGate(gate_name.clone()))?;

        let qubits = gc
            .qubits()
            .iter()
            .map(|texpr| extract_qubit_index(texpr, &qubit_base))
            .collect::<Result<Vec<_>, _>>()?;

        let params = match gc.params() {
            None => vec![],
            Some(plist) => plist
                .iter()
                .map(|texpr| eval_texpr(texpr, &symbol_table))
                .collect::<Result<Vec<_>, _>>()?,
        };

        let mut op = Operation::new(kind, qubits, params);
        op.source_framework = Some("openqasm3".to_string());
        op.original_gate_name = Some(gate_name);
        circuit.push(op);
    }

    circuit.validate().map_err(|e| ParseError::SyntaxError(format!("{e:?}")))?;
    Ok(circuit)
}

/// Extracts a flat qubit index from a gate operand expression.
fn extract_qubit_index(
    texpr: &oq3_semantics::asg::TExpr,
    qubit_base: &HashMap<SymbolId, usize>,
) -> Result<usize, ParseError> {
    let go = match texpr.expression() {
        Expr::GateOperand(go) => go,
        other => {
            return Err(ParseError::SyntaxError(format!(
                "expected gate operand, got {other:?}"
            )))
        }
    };
    match go {
        asg::GateOperand::Identifier(sym_id_result) => {
            let sym_id = sym_id_result
                .as_ref()
                .map_err(|_| ParseError::UndeclaredQubit("<error>".to_string()))?;
            qubit_base
                .get(sym_id)
                .copied()
                .ok_or_else(|| ParseError::UndeclaredQubit(format!("{sym_id:?}")))
        }
        asg::GateOperand::IndexedIdentifier(ii) => {
            let sym_id = ii
                .identifier()
                .as_ref()
                .map_err(|_| ParseError::UndeclaredQubit("<error>".to_string()))?;
            let base = qubit_base
                .get(sym_id)
                .copied()
                .ok_or_else(|| ParseError::UndeclaredQubit(format!("{sym_id:?}")))?;
            let idx = extract_index_from_operator(ii.indexes())?;
            Ok(base + idx)
        }
        asg::GateOperand::HardwareQubit(hq) => Err(ParseError::SyntaxError(format!(
            "hardware qubit '{}' not supported in circuit translation",
            hq.identifier()
        ))),
    }
}

/// Extracts an integer index from the first index operator of an `IndexedIdentifier`.
fn extract_index_from_operator(
    indexes: &[asg::IndexOperator],
) -> Result<usize, ParseError> {
    let op = indexes
        .first()
        .ok_or_else(|| ParseError::SyntaxError("missing qubit index".to_string()))?;
    match op {
        asg::IndexOperator::ExpressionList(el) => {
            let texpr = el
                .expressions
                .first()
                .ok_or_else(|| ParseError::SyntaxError("empty index list".to_string()))?;
            match texpr.expression() {
                Expr::Literal(Literal::Int(il)) => Ok(*il.value() as usize),
                other => Err(ParseError::SyntaxError(format!(
                    "non-integer qubit index: {other:?}"
                ))),
            }
        }
        asg::IndexOperator::SetExpression(_) => Err(ParseError::SyntaxError(
            "set-expression qubit indexing is not supported".to_string(),
        )),
    }
}

/// Evaluates a typed expression to an `f64` parameter value.
fn eval_texpr(
    texpr: &oq3_semantics::asg::TExpr,
    symbol_table: &SymbolTable,
) -> Result<f64, ParseError> {
    match texpr.expression() {
        Expr::Literal(Literal::Float(fl)) => fl
            .value()
            .parse::<f64>()
            .map_err(|_| ParseError::InvalidParam(fl.value().to_string())),
        Expr::Literal(Literal::Int(il)) => Ok(*il.value() as f64),
        Expr::Identifier(sym_id_result) => {
            let sym_id = sym_id_result
                .as_ref()
                .map_err(|_| ParseError::InvalidParam("<error sym>".to_string()))?;
            let name = symbol_table[sym_id].name();
            builtin_const(name).ok_or_else(|| ParseError::InvalidParam(name.to_string()))
        }
        Expr::UnaryExpr(ue) => {
            let operand = eval_texpr(ue.operand(), symbol_table)?;
            Ok(match ue.op() {
                UnaryOp::Minus => -operand,
                _ => operand,
            })
        }
        Expr::BinaryExpr(be) => {
            let left = eval_texpr(be.left(), symbol_table)?;
            let right = eval_texpr(be.right(), symbol_table)?;
            match be.op() {
                BinaryOp::ArithOp(ArithOp::Add) => Ok(left + right),
                BinaryOp::ArithOp(ArithOp::Sub) => Ok(left - right),
                BinaryOp::ArithOp(ArithOp::Mul) => Ok(left * right),
                BinaryOp::ArithOp(ArithOp::Div) => Ok(left / right),
                other => Err(ParseError::InvalidParam(format!(
                    "unsupported binary op in parameter: {other:?}"
                ))),
            }
        }
        other => Err(ParseError::InvalidParam(format!(
            "unsupported param expression: {other:?}"
        ))),
    }
}

/// Maps built-in OpenQASM 3 constant names to their `f64` values.
fn builtin_const(name: &str) -> Option<f64> {
    match name {
        "pi" => Some(std::f64::consts::PI),
        "euler" | "e" => Some(std::f64::consts::E),
        "tau" => Some(std::f64::consts::TAU),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cforge_core::GateKind;

    const BELL_QASM3: &str = r#"
OPENQASM 3.0;
include "stdgates.inc";
qubit[2] q;
h q[0];
cx q[0], q[1];
"#;

    #[test]
    fn bell_state_qasm3() {
        let circuit = parse_qasm3(BELL_QASM3).unwrap();
        assert_eq!(circuit.num_qubits(), 2);
        assert_eq!(circuit.operations.len(), 2);
        assert_eq!(circuit.operations[0].kind, GateKind::H);
        assert_eq!(circuit.operations[0].qubits, vec![0]);
        assert_eq!(circuit.operations[1].kind, GateKind::Cx);
        assert_eq!(circuit.operations[1].qubits, vec![0, 1]);
    }

    #[test]
    fn provenance_metadata_is_qasm3() {
        let circuit = parse_qasm3(BELL_QASM3).unwrap();
        for op in &circuit.operations {
            assert_eq!(op.source_framework.as_deref(), Some("openqasm3"));
        }
    }

    #[test]
    fn three_qubit_two_registers_qasm3() {
        let source = r#"
OPENQASM 3.0;
include "stdgates.inc";
qubit[2] q;
qubit[1] a;
h q[0];
cx q[0], q[1];
x a[0];
"#;
        let circuit = parse_qasm3(source).unwrap();
        assert_eq!(circuit.num_qubits(), 3);
        assert_eq!(circuit.operations.len(), 3);
        assert_eq!(circuit.operations[2].kind, GateKind::X);
        assert_eq!(circuit.operations[2].qubits, vec![2]);
    }
}
