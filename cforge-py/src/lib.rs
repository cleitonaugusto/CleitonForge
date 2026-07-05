use std::collections::HashMap;
use std::path::Path;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use cforge_backends::{
    BackendError, DEFAULT_SEED, DensityMatrixBackend, NativeStateVectorBackend, NoisyConfig,
    NoisyStatevectorBackend, Q1tSimBackend, QuantRS2Backend, RoqoqoBackend,
    SimulationBackend, SimulationResult,
};
use cforge_backends::noise::NoiseChannel;
use cforge_core::{Circuit, GateKind, Operation};
use cforge_metrics::compute_stats;
use cforge_parser::{normalize_convention, parse_qasm2, parse_qasm3, RzConvention};

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
        "roqoqo" => Ok(Box::new(RoqoqoBackend)),
        "noisy" => Ok(Box::new(NoisyStatevectorBackend {
            config: NoisyConfig::depolarizing(0.001, 0.01),
        })),
        "q1tsim" => Ok(Box::new(Q1tSimBackend)),
        "density-matrix" | "dm" => Ok(Box::new(DensityMatrixBackend::default())),
        other => Err(PyValueError::new_err(format!(
            "unknown backend '{other}'; available: statevector, quantrs2, roqoqo, noisy, q1tsim, density-matrix"
        ))),
    }
}

fn build_noisy_backend(
    depolarizing_1q: Option<f64>,
    depolarizing_2q: Option<f64>,
    amplitude_damping: Option<f64>,
    readout_error: Option<f64>,
) -> Box<dyn SimulationBackend> {
    let single_qubit = depolarizing_1q.map(|p| NoiseChannel::Depolarizing { p })
        .or_else(|| amplitude_damping.map(|g| NoiseChannel::AmplitudeDamping { gamma: g }));
    let two_qubit = depolarizing_2q.map(|p| NoiseChannel::Depolarizing { p })
        .or_else(|| amplitude_damping.map(|g| NoiseChannel::AmplitudeDamping { gamma: g }));
    Box::new(NoisyStatevectorBackend {
        config: NoisyConfig { single_qubit, two_qubit, readout: readout_error },
    })
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

// ── Helpers ───────────────────────────────────────────────────────────────────

fn convention_from_str(s: &str) -> PyResult<RzConvention> {
    match s {
        "standard" | "ibm" | "qiskit" => Ok(RzConvention::Standard),
        "reversed" | "quantrs2"       => Ok(RzConvention::Reversed),
        other => Err(PyValueError::new_err(format!(
            "unknown convention '{other}'; use 'standard' (IBM/Qiskit) or 'reversed' (quantrs2)"
        ))),
    }
}

// ── Public functions ──────────────────────────────────────────────────────────

/// Simulate a circuit and return measurement results.
///
/// Args:
///     circuit (Circuit): Circuit built with the builder API.
///     shots (int): Measurement samples. ``0`` returns exact statevector only.
///     seed (int): PRNG seed for reproducible counts. Default: ``cforge.DEFAULT_SEED``.
///     backend (str): ``"statevector"`` (default), ``"quantrs2"``, ``"roqoqo"``, or ``"noisy"``.
///     depolarizing_1q (float | None): Single-qubit depolarizing error rate (0–1).
///         If set, activates the noisy backend regardless of ``backend`` kwarg.
///     depolarizing_2q (float | None): Two-qubit depolarizing error rate (0–1).
///     amplitude_damping (float | None): T1 decay rate γ per gate step.
///     readout_error (float | None): Classical bit-flip probability at measurement.
///
/// Returns:
///     RunResult: Simulation result with ``.counts``, ``.probabilities``, ``.statevector``.
///
/// Example (noiseless)::
///
///     c = cforge.Circuit(2)
///     c.h(0); c.cx(0, 1)
///     r = cforge.run(c, shots=1024)
///     print(r.top_states(2))   # [("11", 0.5), ("00", 0.5)]
///
/// Example (NISQ noise)::
///
///     r = cforge.run(c, shots=4096, depolarizing_1q=0.001, depolarizing_2q=0.01,
///                    readout_error=0.02)
///     print(r.counts)   # dominated by "00"/"11", some leakage to "01"/"10"
#[pyfunction]
#[pyo3(signature = (
    circuit,
    shots = 0,
    seed = DEFAULT_SEED,
    backend = "statevector",
    depolarizing_1q = None,
    depolarizing_2q = None,
    amplitude_damping = None,
    readout_error = None,
))]
#[allow(clippy::too_many_arguments)]
fn run(
    circuit: &PyCircuit,
    shots: usize,
    seed: u64,
    backend: &str,
    depolarizing_1q:   Option<f64>,
    depolarizing_2q:   Option<f64>,
    amplitude_damping: Option<f64>,
    readout_error:     Option<f64>,
) -> PyResult<PyRunResult> {
    let b: Box<dyn SimulationBackend> = if depolarizing_1q.is_some()
        || depolarizing_2q.is_some()
        || amplitude_damping.is_some()
        || readout_error.is_some()
    {
        build_noisy_backend(depolarizing_1q, depolarizing_2q, amplitude_damping, readout_error)
    } else {
        backend_for(backend)?
    };
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
///     backend (str): ``"statevector"``, ``"quantrs2"``, ``"roqoqo"``, or ``"noisy"``.
///     depolarizing_1q (float | None): Single-qubit depolarizing error rate.
///     depolarizing_2q (float | None): Two-qubit depolarizing error rate.
///     amplitude_damping (float | None): T1 amplitude damping rate γ per gate step.
///     readout_error (float | None): Classical readout bit-flip probability.
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
#[pyo3(signature = (
    source,
    shots = 0,
    seed = DEFAULT_SEED,
    backend = "statevector",
    depolarizing_1q = None,
    depolarizing_2q = None,
    amplitude_damping = None,
    readout_error = None,
))]
#[allow(clippy::too_many_arguments)]
fn run_qasm(
    source: &str,
    shots: usize,
    seed: u64,
    backend: &str,
    depolarizing_1q:   Option<f64>,
    depolarizing_2q:   Option<f64>,
    amplitude_damping: Option<f64>,
    readout_error:     Option<f64>,
) -> PyResult<PyRunResult> {
    let circuit = parse_from_str(source)?;
    let b: Box<dyn SimulationBackend> = if depolarizing_1q.is_some()
        || depolarizing_2q.is_some()
        || amplitude_damping.is_some()
        || readout_error.is_some()
    {
        build_noisy_backend(depolarizing_1q, depolarizing_2q, amplitude_damping, readout_error)
    } else {
        backend_for(backend)?
    };
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

/// Rewrite a circuit's Rz-family gate angles to match a target convention.
///
/// Returns a **new** :class:`Circuit`; the original is unchanged.
///
/// Conventions
/// -----------
/// ``"standard"`` (aliases: ``"ibm"``, ``"qiskit"``)
///     ``Rz(λ) = [[e^{-iλ/2}, 0], [0, e^{+iλ/2}]]``
///     Used by: CleitonForge native, roqoqo, q1tsim, IBM Quantum hardware.
///
/// ``"reversed"`` (alias: ``"quantrs2"``)
///     ``Rz(λ) = [[e^{+iλ/2}, 0], [0, e^{-iλ/2}]]``
///     Used by: quantrs2-core.
///
/// Gates affected
/// --------------
/// ``rz``, ``p`` (phase), ``crz``, ``cp``, ``u`` (phi/lambda only), ``cu`` (phi/lambda/gamma).
/// ``rx``, ``ry``, ``h``, ``cx`` and all Clifford gates are left unchanged.
///
/// Args:
///     circuit (Circuit): Source circuit.
///     from_convention (str): Convention the circuit was written for. Default: ``"reversed"``.
///     to_convention (str): Target convention. Default: ``"standard"``.
///
/// Returns:
///     Circuit: New normalized circuit.
///
/// Example::
///
///     import cforge, math
///
///     # QAOA circuit written assuming quantrs2 convention
///     c = cforge.Circuit(2)
///     c.h(0); c.h(1)
///     c.cx(0, 1); c.rz(-3 * math.pi / 2, 1); c.cx(0, 1)
///     c.rx(-math.pi / 4, 0); c.rx(-math.pi / 4, 1)
///
///     fixed = cforge.normalize(c)                   # Reversed → Standard
///     r_native  = cforge.run(c,     backend="statevector")
///     r_q2_raw  = cforge.run(c,     backend="quantrs2")
///     r_q2_norm = cforge.run(fixed, backend="quantrs2")
///
///     # fidelity helper
///     def fidelity(a, b):
///         s = sum(complex(ar, ai).conjugate() * complex(br, bi)
///                 for (ar, ai), (br, bi) in zip(a, b))
///         return abs(s) ** 2
///
///     print(fidelity(r_native.statevector, r_q2_raw.statevector))   # ≈ 0.0
///     print(fidelity(r_native.statevector, r_q2_norm.statevector))  # ≈ 1.0
#[pyfunction]
#[pyo3(signature = (circuit, from_convention = "reversed", to_convention = "standard"))]
fn normalize(
    circuit: &PyCircuit,
    from_convention: &str,
    to_convention: &str,
) -> PyResult<PyCircuit> {
    let from = convention_from_str(from_convention)?;
    let to   = convention_from_str(to_convention)?;
    Ok(PyCircuit { inner: normalize_convention(&circuit.inner, from, to) })
}

/// Parse an OpenQASM string and normalize its Rz convention in one step.
///
/// Equivalent to ``cforge.normalize(cforge.Circuit.from_qasm(source), ...)``.
///
/// Args:
///     source (str): OpenQASM 2.0 or 3.0 source.
///     from_convention (str): Convention assumed by the source. Default: ``"reversed"``.
///     to_convention (str): Target convention. Default: ``"standard"``.
///
/// Returns:
///     Circuit: Parsed and normalized circuit.
///
/// Example::
///
///     qaoa_qasm = """
///     OPENQASM 2.0;
///     qreg q[2];
///     h q[0]; h q[1];
///     cx q[0],q[1]; rz(-4.71238898) q[1]; cx q[0],q[1];
///     rx(-0.39269908) q[0]; rx(-0.39269908) q[1];
///     """
///     fixed_circuit = cforge.normalize_qasm(qaoa_qasm)
///     r = cforge.run(fixed_circuit, backend="quantrs2")
#[pyfunction]
#[pyo3(signature = (source, from_convention = "reversed", to_convention = "standard"))]
fn normalize_qasm(
    source: &str,
    from_convention: &str,
    to_convention: &str,
) -> PyResult<PyCircuit> {
    let circuit = parse_from_str(source)?;
    let from = convention_from_str(from_convention)?;
    let to   = convention_from_str(to_convention)?;
    Ok(PyCircuit { inner: normalize_convention(&circuit, from, to) })
}

// ── Module ────────────────────────────────────────────────────────────────────

#[pymodule]
fn cleitonforge(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyCircuit>()?;
    m.add_class::<PyRunResult>()?;
    m.add_function(wrap_pyfunction!(run, m)?)?;
    m.add_function(wrap_pyfunction!(run_qasm, m)?)?;
    m.add_function(wrap_pyfunction!(validate_qasm, m)?)?;
    m.add_function(wrap_pyfunction!(normalize, m)?)?;
    m.add_function(wrap_pyfunction!(normalize_qasm, m)?)?;
    m.add("DEFAULT_SEED", DEFAULT_SEED)?;
    m.add("CONVENTION_STANDARD", "standard")?;
    m.add("CONVENTION_REVERSED", "reversed")?;
    Ok(())
}
