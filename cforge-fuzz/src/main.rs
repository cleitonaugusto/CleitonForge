//! cforge-fuzz CLI — phase-sensitive differential fuzzing campaigns.
//!
//! Examples:
//!   cforge-fuzz --device quantrs2                    # hunt real bugs
//!   cforge-fuzz --device conjugated-all --oracle probability
//!       # theorem demo: finds nothing, provably
//!   cforge-fuzz --device conjugated-all --oracle amplitude
//!       # same device, N1 oracle: 2-gate witness in seconds
//!   cforge-fuzz --device conjugated-rz --oracle probability --zoo-dir bug-zoo
//!       # rediscovers the H·Rz·S·H witness class automatically

use std::path::PathBuf;

use clap::{Parser, ValueEnum};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use cforge_backends::{
    ConjugatedStateVectorBackend, ConjugationScope, NativeStateVectorBackend, QuantRS2Backend,
    RoqoqoBackend, SimulationBackend,
};
use cforge_core::Circuit;
use cforge_fuzz::generator::random_circuit;
use cforge_fuzz::oracle::{distance, OracleLevel};
use cforge_fuzz::shrinker::shrink;
use cforge_fuzz::triage::triage;
use cforge_fuzz::zoo::{to_qasm2, ZooEntry};

#[derive(Clone, Copy, ValueEnum)]
enum DeviceArg {
    /// quantrs2 gate matrices (real target; known Rz sign bug in v0.2.0)
    Quantrs2,
    /// roqoqo gate matrices (real target)
    Roqoqo,
    /// Reference simulator with every gate conjugated (theorem demo)
    ConjugatedAll,
    /// Reference simulator with Rz/Phase-family sign bug (quantrs2 class)
    ConjugatedRz,
}

impl DeviceArg {
    fn backend(&self) -> Box<dyn SimulationBackend> {
        match self {
            DeviceArg::Quantrs2 => Box::new(QuantRS2Backend),
            DeviceArg::Roqoqo => Box::new(RoqoqoBackend),
            DeviceArg::ConjugatedAll => Box::new(ConjugatedStateVectorBackend::new(
                ConjugationScope::AllGates,
            )),
            DeviceArg::ConjugatedRz => Box::new(ConjugatedStateVectorBackend::new(
                ConjugationScope::RotationSignOnly,
            )),
        }
    }
}

#[derive(Clone, Copy, ValueEnum)]
enum OracleArg {
    /// N1 — amplitude equality modulo global phase (phase-sensitive)
    Amplitude,
    /// N2 — probability equality (what sampling benchmarks see)
    Probability,
    /// N3 — per-qubit ⟨Z⟩ equality
    Observable,
}

impl From<OracleArg> for OracleLevel {
    fn from(a: OracleArg) -> Self {
        match a {
            OracleArg::Amplitude => OracleLevel::Amplitude,
            OracleArg::Probability => OracleLevel::Probability,
            OracleArg::Observable => OracleLevel::Observable,
        }
    }
}

#[derive(Parser)]
#[command(
    name = "cforge-fuzz",
    about = "Quantum differential fuzzer — finds gate convention divergences automatically"
)]
struct Cli {
    /// Device under test, compared against the native reference
    #[arg(long, value_enum, default_value = "quantrs2")]
    device: DeviceArg,

    /// Oracle strictness level
    #[arg(long, value_enum, default_value = "amplitude")]
    oracle: OracleArg,

    /// Number of random circuits to generate and test
    #[arg(long, default_value_t = 10_000)]
    iterations: usize,

    /// Maximum qubits per generated circuit
    #[arg(long, default_value_t = 3)]
    max_qubits: usize,

    /// Maximum gates per generated circuit
    #[arg(long, default_value_t = 8)]
    max_depth: usize,

    /// PRNG seed (0 = random)
    #[arg(long, default_value_t = 0)]
    seed: u64,

    /// Divergence tolerance
    #[arg(long, default_value_t = 1e-6)]
    tol: f64,

    /// Stop after this many divergences (0 = never stop early)
    #[arg(long, default_value_t = 1)]
    stop_after: usize,

    /// Directory to write bug-zoo JSON entries (omit to skip)
    #[arg(long)]
    zoo_dir: Option<PathBuf>,
}

fn diverges(
    circuit: &Circuit,
    device: &dyn SimulationBackend,
    level: OracleLevel,
    tol: f64,
) -> Option<f64> {
    let a = NativeStateVectorBackend.run(circuit, 0, 0).ok()?.statevector;
    let b = device.run(circuit, 0, 0).ok()?.statevector;
    let d = distance(&a, &b, level);
    (d > tol).then_some(d)
}

fn main() {
    let cli = Cli::parse();
    let device = cli.device.backend();
    let level: OracleLevel = cli.oracle.into();

    let seed = if cli.seed == 0 {
        rand::thread_rng().gen()
    } else {
        cli.seed
    };

    println!("cforge-fuzz — phase-sensitive differential fuzzing");
    println!("  reference : statevector-native");
    println!("  device    : {}", device.name());
    println!("  oracle    : {}", level.label());
    println!("  circuits  : {} (≤{}q, ≤{} gates)", cli.iterations, cli.max_qubits, cli.max_depth);
    println!("  seed      : 0x{seed:x}");
    println!();

    let started = std::time::Instant::now();
    let mut found = 0usize;
    let mut tried = 0usize;

    for i in 0..cli.iterations {
        tried = i + 1;
        let mut rng = StdRng::seed_from_u64(seed.wrapping_add(i as u64));
        let nq = rng.gen_range(1..=cli.max_qubits);
        let depth = rng.gen_range(1..=cli.max_depth);
        let circuit = random_circuit(&mut rng, nq, depth);

        let Some(dist) = diverges(&circuit, device.as_ref(), level, cli.tol) else {
            continue;
        };
        found += 1;

        let minimal = shrink(&circuit, &|cand| {
            diverges(cand, device.as_ref(), level, cli.tol).is_some()
        });

        let sv_ref = NativeStateVectorBackend
            .run(&minimal, 0, 0)
            .expect("reference must simulate the minimal witness")
            .statevector;
        let sv_dev = device
            .run(&minimal, 0, 0)
            .expect("device must simulate the minimal witness")
            .statevector;
        let verdict = triage(&minimal, &sv_ref, &sv_dev, cli.tol);

        println!("── divergence #{found} (circuit {i}, distance {dist:.3e}) ─────────────");
        println!("  minimal witness ({} gates):", minimal.operations.len());
        for line in to_qasm2(&minimal).lines().skip(3) {
            println!("    {line}");
        }
        println!("  class      : {}", verdict.class.label());
        println!("  N1 distance: {:.3e}", verdict.amplitude_distance);
        println!("  N2 distance: {:.3e}", verdict.probability_distance);
        println!("  visibility : {}", verdict.visibility_label());

        if let Some(dir) = &cli.zoo_dir {
            let entry = ZooEntry {
                id: format!("{}-{}-{:03}", device.name(), level.label(), found),
                reference_backend: "statevector-native",
                device_backend: device.name(),
                seed,
                circuit: &minimal,
                triage: &verdict,
            };
            match entry.save(dir) {
                Ok(path) => println!("  zoo entry  : {}", path.display()),
                Err(e) => eprintln!("  zoo write failed: {e}"),
            }
        }
        println!();

        if cli.stop_after > 0 && found >= cli.stop_after {
            break;
        }
    }

    let secs = started.elapsed().as_secs_f64();
    println!("Result: {found} divergence(s) in {tried} circuits ({secs:.1}s).");
    if found > 0 {
        std::process::exit(1);
    }
}
