//! OpenQASM 2.0 parser wrapping the `qasm` crate.
//!
//! Translates a QASM 2 source string into CleitonForge's canonical IR. Only
//! gate applications and qubit register declarations are translated;
//! classical registers, measurements, barriers, and conditionals are
//! silently skipped (they have no bearing on the unitary part of the
//! circuit that backends simulate).
//!
//! Custom gate definitions (`gate mygate(params) qubits { body }`) are fully
//! supported: they are collected on the first pass and inlined recursively
//! at every call site, with qubit and parameter substitution applied.

use std::collections::HashMap;
use std::path::Path;

use cforge_core::{Circuit, Operation};

use crate::error::ParseError;
use crate::translate::translate_gate;

// ── User-defined gate representation ─────────────────────────────────────────

struct UserGateDef {
    qubit_names: Vec<String>,
    angle_names: Vec<String>,
    body:        Vec<qasm::AstNode>,
}

type UserGateMap = HashMap<String, UserGateDef>;

/// Collects all `gate` declarations from the AST into a lookup map.
fn collect_user_gates(ast: &[qasm::AstNode]) -> UserGateMap {
    let mut map = UserGateMap::new();
    for node in ast {
        if let qasm::AstNode::Gate(name, qubit_names, angle_names, body) = node {
            map.insert(
                name.clone(),
                UserGateDef {
                    qubit_names: qubit_names.clone(),
                    angle_names: angle_names.clone(),
                    body:        body.clone(),
                },
            );
        }
    }
    map
}

// ── Main entry point ──────────────────────────────────────────────────────────

/// Parses an OpenQASM 2.0 source string and returns the canonical circuit.
///
/// `search_dir` is the directory used to resolve `include` statements (pass
/// `std::env::current_dir().unwrap()` when there are none).
pub fn parse_qasm2(source: &str, search_dir: &Path) -> Result<Circuit, ParseError> {
    let processed = qasm::process(source, search_dir);
    let mut tokens = qasm::lex(&processed);
    let ast = qasm::parse(&mut tokens)
        .map_err(|e| ParseError::SyntaxError(format!("{e:?}")))?;

    // First pass: collect qubit register declarations and user-defined gates.
    let mut reg_map: HashMap<String, (usize, usize)> = HashMap::new();
    let mut total_qubits: usize = 0;

    for node in &ast {
        if let qasm::AstNode::QReg(name, size) = node {
            let size = *size as usize;
            reg_map.insert(name.clone(), (total_qubits, size));
            total_qubits += size;
        }
    }

    let user_gates = collect_user_gates(&ast);
    let mut circuit = Circuit::new(total_qubits);

    // Second pass: translate gate applications.
    for node in &ast {
        let qasm::AstNode::ApplyGate(gate_name, args, raw_params) = node else {
            continue;
        };

        let raw_params: Vec<f64> = raw_params
            .iter()
            .map(|s| eval_param_with_subst(s, &HashMap::new()))
            .collect::<Result<_, _>>()?;

        let arg_lists: Vec<Vec<usize>> = args
            .iter()
            .map(|arg| resolve_arg_from_reg(arg, &reg_map))
            .collect::<Result<_, _>>()?;

        let broadcast_n = arg_lists.iter().map(|v| v.len()).max().unwrap_or(1);

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

            apply_gate_call(
                gate_name,
                &qubits,
                &raw_params,
                &user_gates,
                &mut circuit,
                0,
            )?;
        }
    }

    circuit.validate().map_err(|e| ParseError::SyntaxError(format!("{e:?}")))?;
    Ok(circuit)
}

// ── Gate application (built-in + user-defined) ────────────────────────────────

const MAX_INLINE_DEPTH: usize = 32;

/// Applies one gate call to the circuit, either as a built-in operation or by
/// recursively inlining a user-defined gate definition.
fn apply_gate_call(
    gate_name: &str,
    qubits:    &[usize],
    params:    &[f64],
    user_gates: &UserGateMap,
    circuit:   &mut Circuit,
    depth:     usize,
) -> Result<(), ParseError> {
    if depth > MAX_INLINE_DEPTH {
        return Err(ParseError::SyntaxError(format!(
            "gate inline depth limit ({MAX_INLINE_DEPTH}) exceeded at '{gate_name}'"
        )));
    }

    if let Some((kind, final_params)) = translate_gate(gate_name, params.to_vec()) {
        let mut op = Operation::new(kind, qubits.to_vec(), final_params);
        op.source_framework = Some("openqasm2".to_string());
        op.original_gate_name = Some(gate_name.to_string());
        circuit.push(op);
        return Ok(());
    }

    if let Some(def) = user_gates.get(gate_name) {
        if qubits.len() != def.qubit_names.len() {
            return Err(ParseError::SyntaxError(format!(
                "gate '{gate_name}' expects {} qubits, got {}",
                def.qubit_names.len(), qubits.len()
            )));
        }
        if params.len() != def.angle_names.len() {
            return Err(ParseError::SyntaxError(format!(
                "gate '{gate_name}' expects {} params, got {}",
                def.angle_names.len(), params.len()
            )));
        }
        return inline_user_gate(def, qubits, params, user_gates, circuit, depth + 1);
    }

    Err(ParseError::UnknownGate(gate_name.to_string()))
}

/// Recursively inlines a user-defined gate by substituting formal qubit/param
/// names with the actual values from the call site.
fn inline_user_gate(
    def:        &UserGateDef,
    act_qubits: &[usize],
    act_params: &[f64],
    user_gates: &UserGateMap,
    circuit:    &mut Circuit,
    depth:      usize,
) -> Result<(), ParseError> {
    let qubit_map: HashMap<&str, usize> = def.qubit_names.iter()
        .zip(act_qubits.iter().copied())
        .map(|(n, q)| (n.as_str(), q))
        .collect();

    let param_map: HashMap<&str, f64> = def.angle_names.iter()
        .zip(act_params.iter().copied())
        .map(|(n, v)| (n.as_str(), v))
        .collect();

    for node in &def.body {
        let qasm::AstNode::ApplyGate(name, args, raw_params) = node else {
            continue;
        };

        let resolved_qubits: Vec<usize> = args.iter()
            .map(|arg| resolve_arg_from_qubit_map(arg, &qubit_map))
            .collect::<Result<_, _>>()?;

        let evaluated_params: Vec<f64> = raw_params.iter()
            .map(|expr| eval_param_with_subst(expr, &param_map))
            .collect::<Result<_, _>>()?;

        apply_gate_call(name, &resolved_qubits, &evaluated_params, user_gates, circuit, depth)?;
    }

    Ok(())
}

// ── Argument resolution ───────────────────────────────────────────────────────

/// Resolves a gate argument to flat qubit indices using the **register map**
/// (for top-level `ApplyGate` nodes that reference actual qubit registers).
fn resolve_arg_from_reg(
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

/// Resolves a gate argument to a single flat qubit index using the **qubit map**
/// (for gate body nodes that reference formal qubit parameter names).
fn resolve_arg_from_qubit_map(
    arg: &qasm::Argument,
    qubit_map: &HashMap<&str, usize>,
) -> Result<usize, ParseError> {
    let name = match arg {
        qasm::Argument::Qubit(n, _) | qasm::Argument::Register(n) => n.as_str(),
    };
    qubit_map.get(name).copied().ok_or_else(|| ParseError::UndeclaredQubit(name.to_string()))
}

// ── Parameter evaluation ──────────────────────────────────────────────────────

/// Evaluates a QASM 2 parameter expression string to an `f64`.
///
/// `param_map` substitutes formal parameter names with their concrete values
/// (used when evaluating expressions inside user-defined gate bodies).
fn eval_param_with_subst(
    expr:      &str,
    param_map: &HashMap<&str, f64>,
) -> Result<f64, ParseError> {
    let s = expr.trim();
    if let Ok(v) = s.parse::<f64>() {
        return Ok(v);
    }
    // Substitute pi and then formal param names, longest first to avoid
    // partial substring matches (e.g. "theta" before "the").
    let mut substituted = s.replace("pi", &std::f64::consts::PI.to_string());
    let mut sorted_params: Vec<(&&str, &f64)> = param_map.iter().collect();
    sorted_params.sort_by_key(|(k, _)| std::cmp::Reverse(k.len()));
    for (name, val) in sorted_params {
        substituted = substituted.replace(name, &val.to_string());
    }
    eval_simple_expr(&substituted).map_err(|_| ParseError::InvalidParam(expr.to_string()))
}

/// Evaluates simple additive/multiplicative expressions over f64 literals.
fn eval_simple_expr(s: &str) -> Result<f64, ()> {
    // Handle +/-
    for (i, c) in s.char_indices().rev() {
        if (c == '+' || c == '-') && i > 0 {
            let left = eval_simple_expr(s[..i].trim())?;
            let right = eval_simple_expr(s[i + 1..].trim())?;
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

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use cforge_core::GateKind;

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
        let source = "OPENQASM 2.0;\nqreg q[2];\nqreg r[2];\ncx q,r;\n";
        let dir = std::env::current_dir().unwrap();
        let circuit = parse_qasm2(source, &dir).unwrap();
        assert_eq!(circuit.num_qubits(), 4);
        assert_eq!(circuit.operations.len(), 2);
        assert_eq!(circuit.operations[0].qubits, vec![0, 2]);
        assert_eq!(circuit.operations[1].qubits, vec![1, 3]);
    }

    #[test]
    fn u2_gate_rewrites_to_u_qasm2() {
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

    // ── Custom gate definition tests ──────────────────────────────────────────

    #[test]
    fn custom_gate_no_params_inlined() {
        // Define CCZ = H·CCX·H and apply it; should expand to 3 operations.
        let source = r#"
OPENQASM 2.0;
gate ccz a, b, c {
    h c;
    ccx a, b, c;
    h c;
}
qreg q[3];
ccz q[0], q[1], q[2];
"#;
        let dir = std::env::current_dir().unwrap();
        let circuit = parse_qasm2(source, &dir).unwrap();
        assert_eq!(circuit.num_qubits(), 3);
        assert_eq!(circuit.operations.len(), 3);
        assert_eq!(circuit.operations[0].kind, GateKind::H);
        assert_eq!(circuit.operations[0].qubits, vec![2]);
        assert_eq!(circuit.operations[1].kind, GateKind::Ccx);
        assert_eq!(circuit.operations[1].qubits, vec![0, 1, 2]);
        assert_eq!(circuit.operations[2].kind, GateKind::H);
        assert_eq!(circuit.operations[2].qubits, vec![2]);
    }

    #[test]
    fn custom_gate_with_params_inlined() {
        // gate rzz(theta) a, b { cx a,b; rz(theta) b; cx a,b; }
        let source = r#"
OPENQASM 2.0;
gate rzz(theta) a, b {
    cx a, b;
    rz(theta) b;
    cx a, b;
}
qreg q[2];
rzz(pi/4) q[0], q[1];
"#;
        let dir = std::env::current_dir().unwrap();
        let circuit = parse_qasm2(source, &dir).unwrap();
        assert_eq!(circuit.operations.len(), 3);
        assert_eq!(circuit.operations[0].kind, GateKind::Cx);
        assert_eq!(circuit.operations[1].kind, GateKind::Rz);
        assert!((circuit.operations[1].params[0] - std::f64::consts::FRAC_PI_4).abs() < 1e-10);
        assert_eq!(circuit.operations[2].kind, GateKind::Cx);
    }

    #[test]
    fn nested_custom_gate_inlined() {
        // gate ccz uses another custom gate (cx_pair) in its body.
        let source = r#"
OPENQASM 2.0;
gate my_cx a, b {
    cx a, b;
}
gate my_bell a, b {
    h a;
    my_cx a, b;
}
qreg q[2];
my_bell q[0], q[1];
"#;
        let dir = std::env::current_dir().unwrap();
        let circuit = parse_qasm2(source, &dir).unwrap();
        assert_eq!(circuit.operations.len(), 2);
        assert_eq!(circuit.operations[0].kind, GateKind::H);
        assert_eq!(circuit.operations[1].kind, GateKind::Cx);
    }

    #[test]
    fn custom_gate_called_twice_independent() {
        let source = r#"
OPENQASM 2.0;
gate ccz a, b, c { h c; ccx a,b,c; h c; }
qreg q[3];
qreg r[3];
ccz q[0], q[1], q[2];
ccz r[0], r[1], r[2];
"#;
        let dir = std::env::current_dir().unwrap();
        let circuit = parse_qasm2(source, &dir).unwrap();
        assert_eq!(circuit.num_qubits(), 6);
        assert_eq!(circuit.operations.len(), 6);
        // Second CCZ maps to qubits 3,4,5
        assert_eq!(circuit.operations[3].qubits, vec![5]); // H on r[2]
        assert_eq!(circuit.operations[4].qubits, vec![3, 4, 5]); // CCX r[0..2]
    }
}
