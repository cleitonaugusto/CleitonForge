pub mod ibm_profile;

use std::path::{Path, PathBuf};

use clap::{Parser, Subcommand, ValueEnum};
use comfy_table::{Table, presets::UTF8_FULL};

use cforge_backends::{DEFAULT_SEED, NativeStateVectorBackend, QuantRS2Backend, SimulationBackend};
use cforge_core::MetricsResult;
use cforge_metrics::{compute_stats, measure};
use cforge_parser::{parse_qasm2, parse_qasm3};

#[derive(Parser)]
#[command(
    name = "cforge",
    about = "CleitonForge — neutral benchmarking layer for quantum simulators",
    version
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Clone, ValueEnum)]
enum OutputFormat {
    Table,
    Json,
}

#[derive(Subcommand)]
enum Commands {
    /// Run a circuit on one or more backends and compare metrics.
    Run {
        /// Path to a .qasm file (OpenQASM 2 or 3)
        #[arg(long)]
        circuit: PathBuf,

        /// Comma-separated list of backends: statevector, quantrs2
        #[arg(long, default_value = "statevector,quantrs2")]
        backends: String,

        /// Number of measurement shots (0 = statevector only)
        #[arg(long, default_value_t = 0)]
        shots: usize,

        /// Output format
        #[arg(long, value_enum, default_value_t = OutputFormat::Table)]
        format: OutputFormat,

        /// Seed for the shot-sampling PRNG (default: 0xdeadbeef_cafebabe).
        /// Setting this to a fixed value makes shot counts reproducible across
        /// runs; changing it lets you measure statistical variance.
        #[arg(long, default_value_t = DEFAULT_SEED)]
        seed: u64,
    },

    /// Parse a circuit and show its statistics without running simulation.
    Validate {
        /// Path to a .qasm file (OpenQASM 2 or 3)
        #[arg(long)]
        circuit: PathBuf,
    },
}

fn main() {
    let cli = Cli::parse();

    match cli.command {
        Commands::Run { circuit, backends, shots, format, seed } => {
            cmd_run(&circuit, &backends, shots, format, seed);
        }
        Commands::Validate { circuit } => cmd_validate(&circuit),
    }
}

// ── cforge run ───────────────────────────────────────────────────────────────

fn cmd_run(path: &PathBuf, backends_str: &str, shots: usize, format: OutputFormat, seed: u64) {
    let source = std::fs::read_to_string(path).unwrap_or_else(|e| {
        eprintln!("error: cannot read {:?}: {e}", path);
        std::process::exit(1);
    });

    let circuit = load_circuit(&source, path);
    let stats = compute_stats(&circuit);

    let selected: Vec<Box<dyn SimulationBackend>> = backends_str
        .split(',')
        .map(|name| -> Box<dyn SimulationBackend> {
            match name.trim() {
                "statevector" | "native" => Box::new(NativeStateVectorBackend),
                "quantrs2" => Box::new(QuantRS2Backend),
                other => {
                    eprintln!("error: unknown backend '{other}'. Available: statevector, quantrs2");
                    std::process::exit(1);
                }
            }
        })
        .collect();

    // Native statevector used as fidelity reference (shots=0 → seed irrelevant).
    let ref_result = NativeStateVectorBackend.run(&circuit, 0, DEFAULT_SEED).ok();
    let reference = ref_result.as_ref().map(|r| r.statevector.as_slice());

    let results: Vec<Result<MetricsResult, _>> = selected
        .iter()
        .map(|b| measure(b.as_ref(), &circuit, shots, seed, reference))
        .collect();

    match format {
        OutputFormat::Table => print_table(&circuit, &stats, shots, seed, &selected, &results),
        OutputFormat::Json => print_json(path, &circuit, &stats, shots, seed, &selected, &results),
    }
}

fn print_table(
    circuit: &cforge_core::Circuit,
    stats: &cforge_metrics::CircuitStats,
    shots: usize,
    seed: u64,
    backends: &[Box<dyn SimulationBackend>],
    results: &[Result<MetricsResult, cforge_backends::BackendError>],
) {
    println!(
        "Circuit: {} qubits  |  {} gates  |  depth {}{}",
        circuit.num_qubits(),
        stats.gate_count,
        stats.depth,
        if shots > 0 { format!("  |  seed 0x{seed:x}") } else { String::new() },
    );

    let mut table = Table::new();
    table.load_preset(UTF8_FULL);
    table.set_header(["Backend", "Time (ms)", "Memory", "Depth", "Gates", "Fidelity", "Shots"]);

    for (backend, result) in backends.iter().zip(results) {
        match result {
            Ok(m) => {
                let fidelity_str = m.fidelity
                    .map(|f| format!("{:.6}", f))
                    .unwrap_or_else(|| "—".to_string());
                let shots_str = if shots > 0 { shots.to_string() } else { "—".to_string() };
                let mem_str = m.memory_bytes.map(format_bytes).unwrap_or_else(|| "—".to_string());
                table.add_row([
                    m.backend_name.as_str(),
                    &format!("{:.3}", m.execution_time_ms),
                    &mem_str,
                    &m.depth.to_string(),
                    &m.gate_count.to_string(),
                    &fidelity_str,
                    &shots_str,
                ]);
            }
            Err(e) => {
                table.add_row([backend.name(), "ERROR", "—", "—", "—", &e.to_string(), "—"]);
            }
        }
    }

    println!("{table}");
}

fn print_json(
    path: &Path,
    circuit: &cforge_core::Circuit,
    stats: &cforge_metrics::CircuitStats,
    shots: usize,
    seed: u64,
    backends: &[Box<dyn SimulationBackend>],
    results: &[Result<MetricsResult, cforge_backends::BackendError>],
) {
    let backend_results: Vec<serde_json::Value> = backends
        .iter()
        .zip(results)
        .map(|(backend, result)| match result {
            Ok(m) => serde_json::json!({
                "backend":      m.backend_name,
                "time_ms":      m.execution_time_ms,
                "memory_bytes": m.memory_bytes,
                "depth":        m.depth,
                "gates":        m.gate_count,
                "fidelity":     m.fidelity,
                "shots":        if shots > 0 { Some(shots) } else { None::<usize> },
                "error":        null,
            }),
            Err(e) => serde_json::json!({
                "backend": backend.name(),
                "error":   e.to_string(),
            }),
        })
        .collect();

    let output = serde_json::json!({
        "circuit": {
            "file":   path.display().to_string(),
            "qubits": circuit.num_qubits(),
            "gates":  stats.gate_count,
            "depth":  stats.depth,
        },
        "shots": shots,
        "seed":  seed,
        "results": backend_results,
    });

    println!("{}", serde_json::to_string_pretty(&output).unwrap());
}

// ── cforge validate ──────────────────────────────────────────────────────────

fn cmd_validate(path: &PathBuf) {
    let source = std::fs::read_to_string(path).unwrap_or_else(|e| {
        eprintln!("error: cannot read {:?}: {e}", path);
        std::process::exit(1);
    });

    let circuit = load_circuit(&source, path);
    let stats = compute_stats(&circuit);

    println!("File   : {}", path.display());
    println!("Qubits : {}", circuit.num_qubits());
    println!("Gates  : {}", stats.gate_count);
    println!("Depth  : {}", stats.depth);

    if !stats.gate_counts_by_kind.is_empty() {
        let mut kinds: Vec<_> = stats.gate_counts_by_kind.iter().collect();
        kinds.sort_by_key(|(k, _)| k.qasm_name());
        println!("By gate:");
        for (kind, count) in kinds {
            println!("  {:8} {}", kind.qasm_name(), count);
        }
    }

    if circuit.validate().is_ok() {
        println!("Status : OK");
    } else {
        println!("Status : INVALID");
        std::process::exit(1);
    }
}

// ── shared helpers ───────────────────────────────────────────────────────────

fn format_bytes(bytes: u64) -> String {
    const GB: u64 = 1 << 30;
    const MB: u64 = 1 << 20;
    const KB: u64 = 1 << 10;
    match bytes {
        b if b >= GB => format!("{:.1} GB", b as f64 / GB as f64),
        b if b >= MB => format!("{:.1} MB", b as f64 / MB as f64),
        b if b >= KB => format!("{:.0} KB", b as f64 / KB as f64),
        b => format!("{} B", b),
    }
}

fn load_circuit(source: &str, path: &PathBuf) -> cforge_core::Circuit {
    let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("");
    let is_qasm3 = source.trim_start().starts_with("OPENQASM 3");

    let result = if is_qasm3 {
        parse_qasm3(source)
    } else {
        let dir = path.parent().unwrap_or(std::path::Path::new("."));
        parse_qasm2(source, dir)
    };

    result.unwrap_or_else(|e| {
        eprintln!("error: failed to parse {:?} ({}): {e}", path, ext);
        std::process::exit(1);
    })
}
