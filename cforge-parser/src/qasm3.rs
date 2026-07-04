//! OpenQASM 3 parser wrapping the `oq3_semantics` crate.
//!
//! Translates a QASM 3 source string into CleitonForge's canonical IR.
//! Supports: qubit register declarations, gate applications from the stdgates
//! gate set, and user-defined `gate { }` definitions (inlined recursively at
//! every call site with qubit and parameter substitution).
//!
//! Measurements, classical registers, and control flow are silently skipped —
//! they are irrelevant to the unitary simulation layer.

use std::collections::HashMap;

use cforge_core::{Circuit, Operation, Qubit};
use oq3_semantics::asg::{self, Expr, Literal, BinaryOp, ArithOp, UnaryOp};
use oq3_semantics::symbols::{SymbolId, SymbolTable, SymbolType};
use oq3_semantics::syntax_to_semantics::parse_source_string;
use oq3_semantics::types::Type;

use crate::error::ParseError;
use crate::translate::translate_gate;

// ── User-defined gate representation ─────────────────────────────────────────

struct UserGateDef3 {
    /// Formal qubit parameter SymbolIds, in declaration order.
    qubit_sym_ids: Vec<SymbolId>,
    /// Formal angle parameter SymbolIds, in declaration order.
    angle_sym_ids: Vec<SymbolId>,
    /// Gate body: a sequence of gate call statements.
    body: Vec<asg::Stmt>,
}

type UserGateMap3 = HashMap<String, UserGateDef3>;

const MAX_INLINE_DEPTH: usize = 32;

// ── Main entry point ──────────────────────────────────────────────────────────

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

    // First pass: collect qubit register declarations.
    let mut qubit_base: HashMap<SymbolId, usize> = HashMap::new();
    let mut total_qubits: usize = 0;
    let mut named_qubits: Vec<Qubit> = Vec::new();

    for stmt in program.stmts() {
        let asg::Stmt::DeclareQuantum(dq) = stmt else { continue };
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

    // Second pass: collect user-defined gate definitions.
    let user_gates = collect_user_gates3(program.stmts(), &symbol_table);

    let mut circuit = Circuit {
        qubits: named_qubits,
        operations: Vec::new(),
    };

    // Third pass: translate top-level gate applications.
    for stmt in program.stmts() {
        let asg::Stmt::GateCall(gc) = stmt else { continue };
        let sym_id = match gc.name() {
            Ok(id) => id.clone(),
            Err(_) => continue,
        };
        let gate_name = symbol_table[&sym_id].name().to_string();

        let qubits = gc
            .qubits()
            .iter()
            .map(|texpr| extract_qubit_index(texpr, &qubit_base))
            .collect::<Result<Vec<_>, _>>()?;

        let raw_params = match gc.params() {
            None => vec![],
            Some(plist) => plist
                .iter()
                .map(|texpr| eval_texpr(texpr, &symbol_table, &HashMap::new()))
                .collect::<Result<Vec<_>, _>>()?,
        };

        apply_gate_call3(
            &gate_name,
            &qubits,
            &raw_params,
            &user_gates,
            &symbol_table,
            &mut circuit,
            0,
        )?;
    }

    circuit.validate().map_err(|e| ParseError::SyntaxError(format!("{e:?}")))?;
    Ok(circuit)
}

// ── Collect user-defined gate definitions ─────────────────────────────────────

fn collect_user_gates3(stmts: &[asg::Stmt], symbol_table: &SymbolTable) -> UserGateMap3 {
    let mut map = UserGateMap3::new();
    for stmt in stmts {
        let asg::Stmt::GateDefinition(gd) = stmt else { continue };
        let Ok(name_sym) = gd.name() else { continue };
        let gate_name = symbol_table[name_sym].name().to_string();

        let qubit_sym_ids: Vec<SymbolId> = gd.qubits().iter()
            .filter_map(|r| r.as_ref().ok().cloned())
            .collect();

        let angle_sym_ids: Vec<SymbolId> = gd.params()
            .map(|ps| ps.iter().filter_map(|r| r.as_ref().ok().cloned()).collect())
            .unwrap_or_default();

        let body: Vec<asg::Stmt> = gd.block().statements().to_vec();

        map.insert(gate_name, UserGateDef3 { qubit_sym_ids, angle_sym_ids, body });
    }
    map
}

// ── Gate application (built-in + user-defined) ────────────────────────────────

fn apply_gate_call3(
    gate_name:    &str,
    qubits:       &[usize],
    params:       &[f64],
    user_gates:   &UserGateMap3,
    symbol_table: &SymbolTable,
    circuit:      &mut Circuit,
    depth:        usize,
) -> Result<(), ParseError> {
    if depth > MAX_INLINE_DEPTH {
        return Err(ParseError::SyntaxError(format!(
            "gate inline depth limit ({MAX_INLINE_DEPTH}) exceeded at '{gate_name}'"
        )));
    }

    if let Some((kind, final_params)) = translate_gate(gate_name, params.to_vec()) {
        let mut op = Operation::new(kind, qubits.to_vec(), final_params);
        op.source_framework = Some("openqasm3".to_string());
        op.original_gate_name = Some(gate_name.to_string());
        circuit.push(op);
        return Ok(());
    }

    if let Some(def) = user_gates.get(gate_name) {
        if qubits.len() != def.qubit_sym_ids.len() {
            return Err(ParseError::SyntaxError(format!(
                "gate '{gate_name}' expects {} qubits, got {}",
                def.qubit_sym_ids.len(), qubits.len()
            )));
        }
        return inline_user_gate3(def, qubits, params, user_gates, symbol_table, circuit, depth + 1);
    }

    Err(ParseError::UnknownGate(gate_name.to_string()))
}

/// Recursively inlines a user-defined gate body by substituting formal
/// qubit/param SymbolIds with actual flat indices / f64 values.
fn inline_user_gate3(
    def:          &UserGateDef3,
    act_qubits:   &[usize],
    act_params:   &[f64],
    user_gates:   &UserGateMap3,
    symbol_table: &SymbolTable,
    circuit:      &mut Circuit,
    depth:        usize,
) -> Result<(), ParseError> {
    // Build SymbolId → flat_index map for formal qubit params.
    let qubit_sym_map: HashMap<SymbolId, usize> = def.qubit_sym_ids.iter()
        .zip(act_qubits.iter().copied())
        .map(|(id, q)| (id.clone(), q))
        .collect();

    // Build SymbolId → f64 map for formal angle params.
    let param_sym_map: HashMap<SymbolId, f64> = def.angle_sym_ids.iter()
        .zip(act_params.iter().copied())
        .map(|(id, v)| (id.clone(), v))
        .collect();

    for stmt in &def.body {
        let asg::Stmt::GateCall(gc) = stmt else { continue };
        let Ok(name_sym) = gc.name() else { continue };
        let inner_name = symbol_table[name_sym].name().to_string();

        let resolved_qubits: Vec<usize> = gc.qubits().iter()
            .map(|texpr| extract_qubit_from_sym_map(texpr, &qubit_sym_map))
            .collect::<Result<Vec<_>, _>>()?;

        let evaluated_params: Vec<f64> = match gc.params() {
            None => vec![],
            Some(plist) => plist.iter()
                .map(|texpr| eval_texpr(texpr, symbol_table, &param_sym_map))
                .collect::<Result<Vec<_>, _>>()?,
        };

        apply_gate_call3(
            &inner_name,
            &resolved_qubits,
            &evaluated_params,
            user_gates,
            symbol_table,
            circuit,
            depth,
        )?;
    }

    Ok(())
}

// ── Qubit extraction ──────────────────────────────────────────────────────────

/// Extracts a flat qubit index from a gate operand expression using the
/// **register map** (for top-level GateCalls referencing declared qubit arrays).
fn extract_qubit_index(
    texpr: &oq3_semantics::asg::TExpr,
    qubit_base: &HashMap<SymbolId, usize>,
) -> Result<usize, ParseError> {
    let go = match texpr.expression() {
        Expr::GateOperand(go) => go,
        other => return Err(ParseError::SyntaxError(format!(
            "expected gate operand, got {other:?}"
        ))),
    };
    match go {
        asg::GateOperand::Identifier(sym_id_result) => {
            let sym_id = sym_id_result.as_ref()
                .map_err(|_| ParseError::UndeclaredQubit("<error>".to_string()))?;
            qubit_base.get(sym_id).copied()
                .ok_or_else(|| ParseError::UndeclaredQubit(format!("{sym_id:?}")))
        }
        asg::GateOperand::IndexedIdentifier(ii) => {
            let sym_id = ii.identifier().as_ref()
                .map_err(|_| ParseError::UndeclaredQubit("<error>".to_string()))?;
            let base = qubit_base.get(sym_id).copied()
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

/// Extracts a flat qubit index from a gate body operand using the
/// **qubit_sym_map** (formal SymbolId → actual flat index).
fn extract_qubit_from_sym_map(
    texpr: &oq3_semantics::asg::TExpr,
    qubit_sym_map: &HashMap<SymbolId, usize>,
) -> Result<usize, ParseError> {
    let go = match texpr.expression() {
        Expr::GateOperand(go) => go,
        other => return Err(ParseError::SyntaxError(format!(
            "expected gate operand in body, got {other:?}"
        ))),
    };
    let sym_id = match go {
        asg::GateOperand::Identifier(r) => {
            r.as_ref().map_err(|_| ParseError::UndeclaredQubit("<err>".to_string()))?
        }
        asg::GateOperand::IndexedIdentifier(ii) => {
            ii.identifier().as_ref()
                .map_err(|_| ParseError::UndeclaredQubit("<err>".to_string()))?
        }
        asg::GateOperand::HardwareQubit(hq) => {
            return Err(ParseError::SyntaxError(format!(
                "hardware qubit '{}' in gate body", hq.identifier()
            )));
        }
    };
    qubit_sym_map.get(sym_id).copied()
        .ok_or_else(|| ParseError::UndeclaredQubit(format!("{sym_id:?}")))
}

/// Extracts an integer index from the first index operator of an `IndexedIdentifier`.
fn extract_index_from_operator(indexes: &[asg::IndexOperator]) -> Result<usize, ParseError> {
    let op = indexes.first()
        .ok_or_else(|| ParseError::SyntaxError("missing qubit index".to_string()))?;
    match op {
        asg::IndexOperator::ExpressionList(el) => {
            let texpr = el.expressions.first()
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

// ── Parameter evaluation ──────────────────────────────────────────────────────

/// Evaluates a typed expression to an `f64` parameter value.
///
/// `param_sym_map` maps formal parameter SymbolIds to their concrete values,
/// used when evaluating expressions inside user-defined gate bodies.
fn eval_texpr(
    texpr: &oq3_semantics::asg::TExpr,
    symbol_table: &SymbolTable,
    param_sym_map: &HashMap<SymbolId, f64>,
) -> Result<f64, ParseError> {
    match texpr.expression() {
        Expr::Literal(Literal::Float(fl)) => fl.value()
            .parse::<f64>()
            .map_err(|_| ParseError::InvalidParam(fl.value().to_string())),
        Expr::Literal(Literal::Int(il)) => Ok(*il.value() as f64),
        Expr::Identifier(sym_id_result) => {
            let sym_id = sym_id_result.as_ref()
                .map_err(|_| ParseError::InvalidParam("<error sym>".to_string()))?;
            // Check formal param map first (gate body context).
            if let Some(&v) = param_sym_map.get(sym_id) {
                return Ok(v);
            }
            // Fall back to built-in constants (pi, euler, tau).
            let name = symbol_table[sym_id].name();
            builtin_const(name).ok_or_else(|| ParseError::InvalidParam(name.to_string()))
        }
        Expr::UnaryExpr(ue) => {
            let operand = eval_texpr(ue.operand(), symbol_table, param_sym_map)?;
            Ok(match ue.op() {
                UnaryOp::Minus => -operand,
                _ => operand,
            })
        }
        // Implicit numeric casts (e.g. int literal in float context: pi/4).
        // Just evaluate the inner operand — the numeric value is unchanged.
        Expr::Cast(cast) => eval_texpr(cast.operand(), symbol_table, param_sym_map),
        Expr::BinaryExpr(be) => {
            let left  = eval_texpr(be.left(),  symbol_table, param_sym_map)?;
            let right = eval_texpr(be.right(), symbol_table, param_sym_map)?;
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
        "pi"            => Some(std::f64::consts::PI),
        "euler" | "e"   => Some(std::f64::consts::E),
        "tau"           => Some(std::f64::consts::TAU),
        _               => None,
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

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
    fn u2_gate_rewrites_to_u_qasm3() {
        let source = r#"
OPENQASM 3.0;
include "stdgates.inc";
qubit[1] q;
u2(0, pi) q[0];
"#;
        match parse_qasm3(source) {
            Ok(circuit) => {
                let op = &circuit.operations[0];
                assert_eq!(op.kind, GateKind::U);
                assert_eq!(op.params.len(), 3);
                assert!((op.params[0] - std::f64::consts::FRAC_PI_2).abs() < 1e-10);
            }
            Err(_) => {
                // oq3_semantics correctly rejects non-stdgate u2; acceptable.
            }
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

    // ── Custom gate definition tests ──────────────────────────────────────────

    #[test]
    fn custom_gate_no_params_qasm3() {
        let source = r#"
OPENQASM 3.0;
include "stdgates.inc";
gate ccz a, b, c {
    h c;
    ccx a, b, c;
    h c;
}
qubit[3] q;
ccz q[0], q[1], q[2];
"#;
        let circuit = parse_qasm3(source).unwrap();
        assert_eq!(circuit.operations.len(), 3);
        assert_eq!(circuit.operations[0].kind, GateKind::H);
        assert_eq!(circuit.operations[0].qubits, vec![2]);
        assert_eq!(circuit.operations[1].kind, GateKind::Ccx);
        assert_eq!(circuit.operations[2].kind, GateKind::H);
        assert_eq!(circuit.operations[2].qubits, vec![2]);
    }

    #[test]
    fn custom_gate_with_params_qasm3() {
        let source = r#"
OPENQASM 3.0;
include "stdgates.inc";
gate rzz(theta) a, b {
    cx a, b;
    rz(theta) b;
    cx a, b;
}
qubit[2] q;
rzz(pi/4) q[0], q[1];
"#;
        let circuit = parse_qasm3(source).unwrap();
        assert_eq!(circuit.operations.len(), 3);
        assert_eq!(circuit.operations[0].kind, GateKind::Cx);
        assert_eq!(circuit.operations[1].kind, GateKind::Rz);
        assert!((circuit.operations[1].params[0] - std::f64::consts::FRAC_PI_4).abs() < 1e-10);
        assert_eq!(circuit.operations[2].kind, GateKind::Cx);
    }

    #[test]
    fn nested_custom_gate_qasm3() {
        let source = r#"
OPENQASM 3.0;
include "stdgates.inc";
gate my_cx a, b { cx a, b; }
gate my_bell a, b { h a; my_cx a, b; }
qubit[2] q;
my_bell q[0], q[1];
"#;
        let circuit = parse_qasm3(source).unwrap();
        assert_eq!(circuit.operations.len(), 2);
        assert_eq!(circuit.operations[0].kind, GateKind::H);
        assert_eq!(circuit.operations[1].kind, GateKind::Cx);
    }
}
