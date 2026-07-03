use std::collections::HashMap;
use std::path::Path;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use cforge_backends::{
    BackendError, DEFAULT_SEED, NativeStateVectorBackend, QuantRS2Backend, SimulationBackend,
    SimulationResult,
};
use cforge_core::{Circuit, GateKind, Operation};
use cforge_metrics::compute_stats;
use cforge_parser::{parse_qasm2, parse_qasm3};

// ── PyCircuit ─────────────────────────────────────────────────────────────────

/// A quantum circuit.
///
/// Build gate-by-gate then pass to :func:`run` or :func:`run_qasm`.
///
/// Example::
///
///     import cforge
///     c = cforge.Circuit(2)
///     c.h(0)
///     c.cx(0, 1)
///     result = cforge.run(c, shots=1024)
///     print(result.counts)   # {"00": 505, "11": 519}
#[pyclass(name = "Circuit")]
pub struct PyCircuit {
    inner: Circuit,
}

impl PyCircuit {
    fn push1(&mut self, kind: GateKind, q: usize) {
        self.inner.push(Operation::new(kind, vec![q], vec![]));
    }
    fn push1p(&mut self, kind: GateKind, q: usize, p: Vec<f64>) {
        self.inner.push(Operation::new(kind, vec![q], p));
    }
    fn push2(&mut self, kind: GateKind, c: usize, t: usize) {
        self.inner.push(Operation::new(kind, vec![c, t], vec![]));
    }
    fn push2p(&mut self, kind: GateKind, c: usize, t: usize, p: Vec<f64>) {
        self.inner.push(Operation::new(kind, vec![c, t], p));
    }
}

#[pymethods]
impl PyCircuit {
    #[new]
    fn new(num_qubits: usize) -> Self {
        Self { inner: Circuit::new(num_qubits) }
    }

    // ── Single-qubit, no params ──────────────────────────────────────────────
    fn id(&mut self, qubit: usize)   { self.push1(GateKind::Id,   qubit); }
    fn x(&mut self, qubit: usize)    { self.push1(GateKind::X,    qubit); }
    fn y(&mut self, qubit: usize)    { self.push1(GateKind::Y,    qubit); }
    fn z(&mut self, qubit: usize)    { self.push1(GateKind::Z,    qubit); }
    fn h(&mut self, qubit: usize)    { self.push1(GateKind::H,    qubit); }
    fn s(&mut self, qubit: usize)    { self.push1(GateKind::S,    qubit); }
    fn sdg(&mut self, qubit: usize)  { self.push1(GateKind::Sdg,  qubit); }
    fn t(&mut self, qubit: usize)    { self.push1(GateKind::T,    qubit); }
    fn tdg(&mut self, qubit: usize)  { self.push1(GateKind::Tdg,  qubit); }
    fn sx(&mut self, qubit: usize)   { self.push1(GateKind::Sx,   qubit); }
    fn sxdg(&mut self, qubit: usize) { self.push1(GateKind::Sxdg, qubit); }

    // ── Single-qubit, parametric ─────────────────────────────────────────────
    fn rx(&mut self, theta: f64, qubit: usize) { self.push1p(GateKind::Rx,    qubit, vec![theta]); }
    fn ry(&mut self, theta: f64, qubit: usize) { self.push1p(GateKind::Ry,    qubit, vec![theta]); }
    fn rz(&mut self, theta: f64, qubit: usize) { self.push1p(GateKind::Rz,    qubit, vec![theta]); }
    fn p(&mut self,  theta: f64, qubit: usize) { self.push1p(GateKind::Phase, qubit, vec![theta]); }

    fn u(&mut self, theta: f64, phi: f64, lam: f64, qubit: usize) {
        self.push1p(GateKind::U, qubit, vec![theta, phi, lam]);
    }

    // ── Two-qubit ────────────────────────────────────────────────────────────
    fn cx(&mut self, ctrl: usize, tgt: usize)   { self.push2(GateKind::Cx,   ctrl, tgt); }
    fn cnot(&mut self, ctrl: usize, tgt: usize) { self.cx(ctrl, tgt); }
    fn cy(&mut self, ctrl: usize, tgt: usize)   { self.push2(GateKind::Cy,   ctrl, tgt); }
    fn cz(&mut self, ctrl: usize, tgt: usize)   { self.push2(GateKind::Cz,   ctrl, tgt); }
    fn ch(&mut self, ctrl: usize, tgt: usize)   { self.push2(GateKind::Ch,   ctrl, tgt); }
    fn csx(&mut self, ctrl: usize, tgt: usize)  { self.push2(GateKind::Csx,  ctrl, tgt); }
    fn swap(&mut self, a: usize, b: usize)      { self.push2(GateKind::Swap, a, b); }

    fn crx(&mut self, theta: f64, ctrl: usize, tgt: usize) { self.push2p(GateKind::Crx, ctrl, tgt, vec![theta]); }
    fn cry(&mut self, theta: f64, ctrl: usize, tgt: usize) { self.push2p(GateKind::Cry, ctrl, tgt, vec![theta]); }
    fn crz(&mut self, theta: f64, ctrl: usize, tgt: usize) { self.push2p(GateKind::Crz, ctrl, tgt, vec![theta]); }
    fn cp(&mut self,  theta: f64, ctrl: usize, tgt: usize) { self.push2p(GateKind::Cp,  ctrl, tgt, vec![theta]); }

    fn cu(&mut self, theta: f64, phi: f64, lam: f64, gamma: f64, ctrl: usize, tgt: usize) {
        self.push2p(GateKind::Cu, ctrl, tgt, vec![theta, phi, lam, gamma]);
    }

    // ── Three-qubit ──────────────────────────────────────────────────────────
    fn ccx(&mut self, c0: usize, c1: usize, tgt: usize) {
        self.inner.push(Operation::new(GateKind::Ccx, vec![c0, c1, tgt], vec![]));
    }
    fn toffoli(&mut self, c0: usize, c1: usize, tgt: usize) { self.ccx(c0, c1, tgt); }

    fn cswap(&mut self, ctrl: usize, a: usize, b: usize) {
        self.inner.push(Operation::new(GateKind::Cswap, vec![ctrl, a, b], vec![]));
    }
    fn fredkin(&mut self, ctrl: usize, a: usize, b: usize) { self.cswap(ctrl, a, b); }

    // ── Properties ───────────────────────────────────────────────────────────
    #[getter]
    fn num_qubits(&self) -> usize { self.inner.num_qubits() }

    #[getter]
    fn gate_count(&self) -> usize { self.inner.operations.len() }

    #[getter]
    fn depth(&self) -> usize { compute_stats(&self.inner).depth }

    fn __len__(&self) -> usize { self.inner.operations.len() }

    fn __repr__(&self) -> String {
        let stats = compute_stats(&self.inner);
        format!(
            "Circuit(qubits={}, gates={}, depth={})",
            self.inner.num_qubits(),
            stats.gate_count,
            stats.depth,
        )
    }
}

// ── PyRunResult ───────────────────────────────────────────────────────────────

/// Result of a quantum circuit simulation.
///
/// Attributes:
///     counts (dict[str, int]): Bitstring measurement counts when shots > 0.
///     probabilities (list[float]): Probability of each basis state (always available).
///
/// Statevector is available as a list of ``(re, im)`` tuples via :attr:`statevector`.
/// Convert to Python complex: ``[complex(re, im) for re, im in result.statevector]``.
#[pyclass(name = "RunResult")]
pub struct PyRunResult {
    #[pyo3(get)]
    counts: HashMap<String, usize>,
    #[pyo3(get)]
    probabilities: Vec<f64>,
    sv: Vec<(f64, f64)>,
    n_qubits: usize,
}

#[pymethods]
impl PyRunResult {
    /// Statevector as a list of ``(real, imag)`` tuples.
    #[getter]
    fn statevector(&self) -> Vec<(f64, f64)> {
        self.sv.clone()
    }

    /// Top ``n`` basis states by probability: ``[(bitstring, probability), ...]``.
    #[pyo3(signature = (n = 5))]
    fn top_states(&self, n: usize) -> Vec<(String, f64)> {
        if self.sv.is_empty() {
            return vec![];
        }
        let width = self.n_qubits;
        let mut v: Vec<(String, f64)> = self
            .probabilities
            .iter()
            .enumerate()
            .map(|(i, &p)| (format!("{:0>width$b}", i, width = width), p))
            .collect();
        v.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        v.truncate(n);
        v
    }

    fn __repr__(&self) -> String {
        if !self.counts.is_empty() {
            let mut sorted: Vec<_> = self.counts.iter().collect();
            sorted.sort_by(|a, b| b.1.cmp(a.1));
            let top = sorted
                .iter()
                .take(3)
                .map(|(k, v)| format!("|{}⟩:{}", k, v))
                .collect::<Vec<_>>()
                .join(", ");
            format!("RunResult(states={}, top=[{}])", self.counts.len(), top)
        } else {
            let peak = self.probabilities.iter().cloned().fold(0.0_f64, f64::max);
            format!("RunResult(statevector, peak_prob={:.4})", peak)
        }
    }
}

// ── Helpers ───────────────────────────────────────────────────────────────────

fn backend_for(name: &str) -> PyResult<Box<dyn SimulationBackend>> {
    match name {
        "statevector" | "native" => Ok(Box::new(NativeStateVectorBackend)),
        "quantrs2" => Ok(Box::new(QuantRS2Backend)),
        other => Err(PyValueError::new_err(format!(
            "unknown backend '{other}'; available: statevector, quantrs2"
        ))),
    }
}

fn into_py_result(
    result: Result<SimulationResult, BackendError>,
    n_qubits: usize,
) -> PyResult<PyRunResult> {
    let r = result.map_err(|e| PyValueError::new_err(e.to_string()))?;
    let probabilities = r.statevector.iter().map(|a| a.norm_sqr()).collect();
    let sv = r.statevector.iter().map(|a| (a.re, a.im)).collect();
    Ok(PyRunResult { counts: r.counts, probabilities, sv, n_qubits })
}

/// Strip ``include`` directives that would require disk access.
/// CleitonForge handles all stdgates.inc gates internally.
fn strip_includes(source: &str) -> String {
    source
        .lines()
        .filter(|l| !l.trim_start().starts_with("include"))
        .collect::<Vec<_>>()
        .join("\n")
}

fn parse_from_str(source: &str) -> PyResult<Circuit> {
    let cleaned = strip_includes(source);
    if cleaned.trim_start().starts_with("OPENQASM 3") {
        parse_qasm3(&cleaned).map_err(|e| PyValueError::new_err(e.to_string()))
    } else {
        parse_qasm2(&cleaned, Path::new(".")).map_err(|e| PyValueError::new_err(e.to_string()))
    }
}

// ── Public functions ──────────────────────────────────────────────────────────

/// Simulate a circuit and return measurement results.
///
/// Args:
///     circuit (Circuit): Circuit built with the builder API.
///     shots (int): Measurement samples. ``0`` returns exact statevector only.
///     seed (int): PRNG seed for reproducible counts. Default: ``cforge.DEFAULT_SEED``.
///     backend (str): ``"statevector"`` (default) or ``"quantrs2"``.
///
/// Returns:
///     RunResult: Simulation result with ``.counts``, ``.probabilities``, ``.statevector``.
///
/// Example::
///
///     c = cforge.Circuit(2)
///     c.h(0); c.cx(0, 1)
///     r = cforge.run(c, shots=1024)
///     print(r.top_states(2))   # [("11", 0.5), ("00", 0.5)]
#[pyfunction]
#[pyo3(signature = (circuit, shots = 0, seed = DEFAULT_SEED, backend = "statevector"))]
fn run(circuit: &PyCircuit, shots: usize, seed: u64, backend: &str) -> PyResult<PyRunResult> {
    let b = backend_for(backend)?;
    let n = circuit.inner.num_qubits();
    into_py_result(b.run(&circuit.inner, shots, seed), n)
}

/// Parse an OpenQASM string and simulate it.
///
/// ``include`` directives are stripped automatically — all gates from
/// ``stdgates.inc`` are handled internally.
///
/// Args:
///     source (str): OpenQASM 2.0 or 3.0 source string.
///     shots (int): Measurement samples. ``0`` = exact statevector only.
///     seed (int): PRNG seed.
///     backend (str): ``"statevector"`` or ``"quantrs2"``.
///
/// Returns:
///     RunResult
///
/// Example::
///
///     qasm = """
///     OPENQASM 2.0;
///     qreg q[2];
///     h q[0];
///     cx q[0],q[1];
///     """
///     r = cforge.run_qasm(qasm, shots=1024)
#[pyfunction]
#[pyo3(signature = (source, shots = 0, seed = DEFAULT_SEED, backend = "statevector"))]
fn run_qasm(source: &str, shots: usize, seed: u64, backend: &str) -> PyResult<PyRunResult> {
    let circuit = parse_from_str(source)?;
    let b = backend_for(backend)?;
    let n = circuit.num_qubits();
    into_py_result(b.run(&circuit, shots, seed), n)
}

/// Parse an OpenQASM string and return circuit statistics without simulating.
///
/// Returns:
///     dict with keys ``qubits``, ``gates``, ``depth``.
#[pyfunction]
fn validate_qasm(source: &str) -> PyResult<HashMap<&'static str, usize>> {
    let circuit = parse_from_str(source)?;
    let stats = compute_stats(&circuit);
    let mut m = HashMap::new();
    m.insert("qubits", circuit.num_qubits());
    m.insert("gates",  stats.gate_count);
    m.insert("depth",  stats.depth);
    Ok(m)
}

// ── Module ────────────────────────────────────────────────────────────────────

#[pymodule]
fn cleitonforge(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyCircuit>()?;
    m.add_class::<PyRunResult>()?;
    m.add_function(wrap_pyfunction!(run, m)?)?;
    m.add_function(wrap_pyfunction!(run_qasm, m)?)?;
    m.add_function(wrap_pyfunction!(validate_qasm, m)?)?;
    m.add("DEFAULT_SEED", DEFAULT_SEED)?;
    Ok(())
}
