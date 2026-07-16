//! The bug zoo: a versioned, public JSON corpus of minimal
//! counterexamples, one file per confirmed divergence.
//!
//! Every entry is fully reproducible: pinned backend names, the exact
//! gate list, the QASM 2 reproducer, the generator seed, and the oracle
//! distances at each level.

use std::fs;
use std::io::Write as _;
use std::path::Path;

use cforge_core::Circuit;

use crate::triage::Triage;

/// A single bug-zoo entry, serialized as JSON by hand to avoid pulling
/// serde into the fuzzer's dependency tree.
pub struct ZooEntry<'a> {
    pub id: String,
    pub reference_backend: &'a str,
    pub device_backend: &'a str,
    pub seed: u64,
    pub circuit: &'a Circuit,
    pub triage: &'a Triage,
}

/// Renders the circuit as an OpenQASM 2 reproducer.
pub fn to_qasm2(circuit: &Circuit) -> String {
    let mut out = String::from("OPENQASM 2.0;\ninclude \"qelib1.inc\";\n");
    out.push_str(&format!("qreg q[{}];\n", circuit.num_qubits()));
    for op in &circuit.operations {
        let name = op.kind.qasm_name();
        let params = if op.params.is_empty() {
            String::new()
        } else {
            format!(
                "({})",
                op.params
                    .iter()
                    .map(|p| format!("{p:.12}"))
                    .collect::<Vec<_>>()
                    .join(",")
            )
        };
        let qubits = op
            .qubits
            .iter()
            .map(|q| format!("q[{q}]"))
            .collect::<Vec<_>>()
            .join(",");
        out.push_str(&format!("{name}{params} {qubits};\n"));
    }
    out
}

fn json_escape(s: &str) -> String {
    s.replace('\\', "\\\\").replace('"', "\\\"").replace('\n', "\\n")
}

impl ZooEntry<'_> {
    pub fn to_json(&self) -> String {
        let gates: Vec<String> = self
            .circuit
            .operations
            .iter()
            .map(|op| {
                format!(
                    r#"{{"gate":"{}","qubits":[{}],"params":[{}]}}"#,
                    op.kind.qasm_name(),
                    op.qubits
                        .iter()
                        .map(|q| q.to_string())
                        .collect::<Vec<_>>()
                        .join(","),
                    op.params
                        .iter()
                        .map(|p| format!("{p:.12}"))
                        .collect::<Vec<_>>()
                        .join(",")
                )
            })
            .collect();

        format!(
            r#"{{
  "id": "{id}",
  "reference_backend": "{reference}",
  "device_backend": "{device}",
  "generator_seed": {seed},
  "num_qubits": {nq},
  "num_gates": {ng},
  "gates": [{gates}],
  "qasm2": "{qasm}",
  "amplitude_distance": {amp:.12e},
  "probability_distance": {prob:.12e},
  "classification": "{class}",
  "sampling_visible": {visible},
  "note": "{note}"
}}
"#,
            id = json_escape(&self.id),
            reference = json_escape(self.reference_backend),
            device = json_escape(self.device_backend),
            seed = self.seed,
            nq = self.circuit.num_qubits(),
            ng = self.circuit.operations.len(),
            gates = gates.join(","),
            qasm = json_escape(&to_qasm2(self.circuit)),
            amp = self.triage.amplitude_distance,
            prob = self.triage.probability_distance,
            class = self.triage.class.label(),
            visible = self.triage.sampling_visible,
            note = if self.triage.sampling_visible {
                "detectable in measurement counts"
            } else {
                "invisible to QV/XEB/mirror/RB/QPE by the conjugation-invariance theorem"
            },
        )
    }

    /// Writes the entry to `<dir>/<id>.json`, creating the directory if
    /// needed.
    pub fn save(&self, dir: &Path) -> std::io::Result<std::path::PathBuf> {
        fs::create_dir_all(dir)?;
        let path = dir.join(format!("{}.json", self.id));
        let mut f = fs::File::create(&path)?;
        f.write_all(self.to_json().as_bytes())?;
        Ok(path)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::triage::{BugClass, Triage};
    use cforge_core::{GateKind, Operation};

    #[test]
    fn qasm2_reproducer_is_well_formed() {
        let mut c = Circuit::new(2);
        c.push(Operation::new(GateKind::H, vec![0], vec![]));
        c.push(Operation::new(GateKind::Rz, vec![0], vec![0.5]));
        c.push(Operation::new(GateKind::Cx, vec![0, 1], vec![]));
        let q = to_qasm2(&c);
        assert!(q.contains("qreg q[2];"));
        assert!(q.contains("rz(0.500000000000) q[0];"));
        assert!(q.contains("cx q[0],q[1];"));
    }

    #[test]
    fn json_round_trip_is_parseable() {
        let mut c = Circuit::new(1);
        c.push(Operation::new(GateKind::H, vec![0], vec![]));
        let t = Triage {
            class: BugClass::RotationConvention,
            amplitude_distance: 1.0,
            probability_distance: 0.0,
            sampling_visible: false,
        };
        let entry = ZooEntry {
            id: "test-001".into(),
            reference_backend: "statevector-native",
            device_backend: "statevector-conjugated",
            seed: 42,
            circuit: &c,
            triage: &t,
        };
        let json = entry.to_json();
        // Minimal structural sanity: balanced braces, key fields present.
        assert_eq!(json.matches('{').count(), json.matches('}').count());
        assert!(json.contains("\"sampling_visible\": false"));
        assert!(json.contains("conjugation-invariance theorem"));
    }
}
