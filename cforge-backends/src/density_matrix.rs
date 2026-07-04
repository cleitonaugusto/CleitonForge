//! Exact noisy simulation via density matrix ρ ∈ C^(2^n × 2^n).
//!
//! Unlike the trajectory-based `NoisyStatevectorBackend` (which is unbiased
//! but has shot noise), this backend tracks the density matrix exactly:
//! a single `run()` call already gives the true probability distribution.
//!
//! ## Representation
//!
//! `rho[i * 2^n + j]` = ρ_{ij}, stored row-major as a flat `Vec<Complex64>`.
//!
//! ## Gate application  ρ → U_k ρ U_k†
//!
//! For qubit k, the unitary U is applied to the row dimension and U† to the
//! column dimension of the 2×2 subspace {i0, i1} × {j0, j1} where only bit k
//! differs. See `apply1_dm` for the exact algorithm.
//!
//! ## Noise (Kraus operators)
//!
//! - Depolarizing: ρ → (1-p)ρ + (p/3)(XρX + YρY + ZρZ)
//! - AmplitudeDamping: ρ → K₀ρK₀† + K₁ρK₁†
//! - BitFlip: ρ → (1-p)ρ + p·XρX
//! - PhaseFlip: ρ → (1-p)ρ + p·ZρZ
//!
//! ## Memory limit
//!
//! At 12 qubits: 4^12 = 16 M Complex64 = 256 MB.
//! At 16 qubits: 4^16 = 4 G entries = 64 GB — not practical.
//! Hard limit: MAX_QUBITS = 12.

use std::collections::HashMap;
use std::f64::consts::{FRAC_1_SQRT_2, PI};

use num_complex::Complex64;

use cforge_core::{Circuit, GateKind};

use crate::noise::{NoisyConfig, NoiseChannel};
use crate::sample::sample_counts;
use crate::trait_def::{BackendError, SimulationBackend, SimulationResult};

const MAX_QUBITS: usize = 12;

/// Exact noisy simulation backend using the density matrix formalism.
///
/// With `config = None` it behaves as a noiseless simulator but tracks ρ
/// rather than |ψ⟩, which is useful for testing and educational purposes.
pub struct DensityMatrixBackend {
    pub config: Option<NoisyConfig>,
}

impl Default for DensityMatrixBackend {
    fn default() -> Self {
        Self { config: None }
    }
}

impl SimulationBackend for DensityMatrixBackend {
    fn name(&self) -> &str {
        "density-matrix"
    }

    fn run(&self, circuit: &Circuit, shots: usize, seed: u64) -> Result<SimulationResult, BackendError> {
        let n = circuit.num_qubits();
        if n > MAX_QUBITS {
            return Err(BackendError(format!(
                "density-matrix backend supports up to {MAX_QUBITS} qubits, circuit has {n}"
            )));
        }

        let big_n = 1usize << n;
        let mut rho = vec![Complex64::new(0.0, 0.0); big_n * big_n];
        rho[0] = Complex64::new(1.0, 0.0); // ρ = |0...0⟩⟨0...0|

        for (op_idx, op) in circuit.operations.iter().enumerate() {
            apply_op(&mut rho, n, op.kind, &op.qubits, &op.params)
                .map_err(|e| BackendError(format!("operation {op_idx}: {e}")))?;

            if let Some(ref cfg) = self.config {
                let is_two_qubit = op.qubits.len() >= 2;
                let noise_ch = if is_two_qubit { cfg.two_qubit.as_ref() } else { cfg.single_qubit.as_ref() };
                if let Some(ch) = noise_ch {
                    for &q in &op.qubits {
                        apply_channel_dm(&mut rho, n, q, ch);
                    }
                }
            }
        }

        // Probabilities from diagonal
        let probs: Vec<f64> = (0..big_n).map(|i| rho[i * big_n + i].re.max(0.0)).collect();
        let probs_sum: f64 = probs.iter().sum();
        let probs: Vec<f64> = probs.into_iter().map(|p| p / probs_sum).collect();

        // Return diagonal as "statevector" (amplitudes = √p, all real)
        let statevector: Vec<Complex64> = probs.iter()
            .map(|&p| Complex64::new(p.sqrt(), 0.0))
            .collect();

        let counts = if shots > 0 {
            // Build a temporary probability-derived pseudostatevector for sample_counts
            sample_counts(&statevector, shots, seed)
        } else {
            HashMap::new()
        };

        Ok(SimulationResult { statevector, counts })
    }
}

// ── Type alias ────────────────────────────────────────────────────────────────

type U2 = [[Complex64; 2]; 2];

// ── Density-matrix gate application ──────────────────────────────────────────

/// Applies ρ → U_k ρ U_k† for qubit k.
fn apply1_dm(rho: &mut [Complex64], n: usize, k: usize, u: &U2) {
    let big_n = 1 << n;
    let bit_k = 1 << k;
    for i0 in 0..big_n {
        if i0 & bit_k != 0 { continue; }
        let i1 = i0 | bit_k;
        for j0 in 0..big_n {
            if j0 & bit_k != 0 { continue; }
            let j1 = j0 | bit_k;

            let m00 = rho[i0 * big_n + j0];
            let m01 = rho[i0 * big_n + j1];
            let m10 = rho[i1 * big_n + j0];
            let m11 = rho[i1 * big_n + j1];

            // T = U * M  (apply U to row dimension)
            let t00 = u[0][0]*m00 + u[0][1]*m10;
            let t01 = u[0][0]*m01 + u[0][1]*m11;
            let t10 = u[1][0]*m00 + u[1][1]*m10;
            let t11 = u[1][0]*m01 + u[1][1]*m11;

            // M' = T * U†  (U†_{ab} = conj(U[b][a]))
            rho[i0 * big_n + j0] = t00*u[0][0].conj() + t01*u[0][1].conj();
            rho[i0 * big_n + j1] = t00*u[1][0].conj() + t01*u[1][1].conj();
            rho[i1 * big_n + j0] = t10*u[0][0].conj() + t11*u[0][1].conj();
            rho[i1 * big_n + j1] = t10*u[1][0].conj() + t11*u[1][1].conj();
        }
    }
}

/// Applies U to `tgt` qubit only within the subspace where `ctrl` = |1⟩.
fn apply_ctrl1_dm(rho: &mut [Complex64], n: usize, ctrl: usize, tgt: usize, u: &U2) {
    let big_n = 1 << n;
    let cb = 1 << ctrl;
    let tb = 1 << tgt;
    for i0 in 0..big_n {
        if (i0 & cb == 0) || (i0 & tb != 0) { continue; }
        let i1 = i0 | tb;
        for j0 in 0..big_n {
            if (j0 & cb == 0) || (j0 & tb != 0) { continue; }
            let j1 = j0 | tb;

            let m00 = rho[i0 * big_n + j0];
            let m01 = rho[i0 * big_n + j1];
            let m10 = rho[i1 * big_n + j0];
            let m11 = rho[i1 * big_n + j1];

            let t00 = u[0][0]*m00 + u[0][1]*m10;
            let t01 = u[0][0]*m01 + u[0][1]*m11;
            let t10 = u[1][0]*m00 + u[1][1]*m10;
            let t11 = u[1][0]*m01 + u[1][1]*m11;

            rho[i0 * big_n + j0] = t00*u[0][0].conj() + t01*u[0][1].conj();
            rho[i0 * big_n + j1] = t00*u[1][0].conj() + t01*u[1][1].conj();
            rho[i1 * big_n + j0] = t10*u[0][0].conj() + t11*u[0][1].conj();
            rho[i1 * big_n + j1] = t10*u[1][0].conj() + t11*u[1][1].conj();
        }
    }
}

/// Applies X to `tgt` only when BOTH `c0` and `c1` are |1⟩ (Toffoli).
fn apply_ccx_dm(rho: &mut [Complex64], n: usize, c0: usize, c1: usize, tgt: usize) {
    let big_n = 1 << n;
    let b0 = 1 << c0; let b1 = 1 << c1; let bt = 1 << tgt;
    for i0 in 0..big_n {
        if (i0 & b0 == 0) || (i0 & b1 == 0) || (i0 & bt != 0) { continue; }
        let i1 = i0 | bt;
        for j0 in 0..big_n {
            if (j0 & b0 == 0) || (j0 & b1 == 0) || (j0 & bt != 0) { continue; }
            let j1 = j0 | bt;
            // X M X† = X M X = swap rows then swap cols
            let m00 = rho[i0 * big_n + j0];
            let m01 = rho[i0 * big_n + j1];
            let m10 = rho[i1 * big_n + j0];
            let m11 = rho[i1 * big_n + j1];
            rho[i0 * big_n + j0] = m11;
            rho[i0 * big_n + j1] = m10;
            rho[i1 * big_n + j0] = m01;
            rho[i1 * big_n + j1] = m00;
        }
    }
}

/// SWAP via CX(a,b)·CX(b,a)·CX(a,b)
fn apply_swap_dm(rho: &mut [Complex64], n: usize, a: usize, b: usize) {
    let x = gate_x();
    apply_ctrl1_dm(rho, n, a, b, &x);
    apply_ctrl1_dm(rho, n, b, a, &x);
    apply_ctrl1_dm(rho, n, a, b, &x);
}

/// CSWAP: apply SWAP(a,b) only when ctrl = |1⟩
fn apply_cswap_dm(rho: &mut [Complex64], n: usize, ctrl: usize, a: usize, b: usize) {
    let big_n = 1 << n;
    let cb = 1 << ctrl; let ba = 1 << a; let bb = 1 << b;
    // Restrict SWAP to the ctrl=1 subspace row-by-row
    // ρ → (|0⟩⟨0|_ctrl ⊗ I + |1⟩⟨1|_ctrl ⊗ SWAP_{ab}) ρ (...)†
    // Build helper: for indices in ctrl=1 subspace, swap a↔b
    for i0 in 0..big_n {
        if (i0 & cb == 0) || (i0 & ba != 0) || (i0 & bb == 0) { continue; }
        let i1 = (i0 | ba) & !bb;
        for j0 in 0..big_n {
            if (j0 & cb == 0) || (j0 & ba != 0) || (j0 & bb == 0) { continue; }
            let j1 = (j0 | ba) & !bb;
            // Swap (i0,j0) ↔ (i0,j1) ↔ (i1,j0) ↔ (i1,j1)
            let a00 = rho[i0 * big_n + j0];
            let a01 = rho[i0 * big_n + j1];
            let a10 = rho[i1 * big_n + j0];
            let a11 = rho[i1 * big_n + j1];
            // SWAP ρ SWAP† = swap both row and col
            rho[i0 * big_n + j0] = a11;
            rho[i0 * big_n + j1] = a10;
            rho[i1 * big_n + j0] = a01;
            rho[i1 * big_n + j1] = a00;
        }
    }
}

// ── Noise channels ────────────────────────────────────────────────────────────

fn apply_channel_dm(rho: &mut [Complex64], n: usize, k: usize, ch: &NoiseChannel) {
    match ch {
        NoiseChannel::Depolarizing { p } => apply_depolarizing_dm(rho, n, k, *p),
        NoiseChannel::AmplitudeDamping { gamma } => apply_amplitude_damping_dm(rho, n, k, *gamma),
        NoiseChannel::BitFlip { p } => {
            let x = gate_x();
            let mut rho_x = rho.to_vec();
            apply1_dm(&mut rho_x, n, k, &x);
            for i in 0..rho.len() {
                rho[i] = (1.0 - *p) * rho[i] + *p * rho_x[i];
            }
        }
        NoiseChannel::PhaseFlip { p } => {
            let z = gate_z();
            let mut rho_z = rho.to_vec();
            apply1_dm(&mut rho_z, n, k, &z);
            for i in 0..rho.len() {
                rho[i] = (1.0 - *p) * rho[i] + *p * rho_z[i];
            }
        }
    }
}

fn apply_depolarizing_dm(rho: &mut [Complex64], n: usize, k: usize, p: f64) {
    let mut rho_x = rho.to_vec();
    let mut rho_y = rho.to_vec();
    let mut rho_z = rho.to_vec();
    apply1_dm(&mut rho_x, n, k, &gate_x());
    apply1_dm(&mut rho_y, n, k, &gate_y());
    apply1_dm(&mut rho_z, n, k, &gate_z());
    for i in 0..rho.len() {
        rho[i] = (1.0 - p) * rho[i] + (p/3.0) * (rho_x[i] + rho_y[i] + rho_z[i]);
    }
}

fn apply_amplitude_damping_dm(rho: &mut [Complex64], n: usize, k: usize, gamma: f64) {
    let s = (1.0 - gamma).sqrt();
    let g = gamma.sqrt();
    let c = Complex64::new;
    // K0 = [[1,0],[0,√(1-γ)]], K1 = [[0,√γ],[0,0]]
    let k0: U2 = [[c(1.0,0.0), c(0.0,0.0)], [c(0.0,0.0), c(s,0.0)]];
    let k1: U2 = [[c(0.0,0.0), c(g,0.0)],  [c(0.0,0.0), c(0.0,0.0)]];
    let mut rho0 = rho.to_vec();
    let mut rho1 = rho.to_vec();
    apply1_dm(&mut rho0, n, k, &k0);
    apply1_dm(&mut rho1, n, k, &k1);
    for i in 0..rho.len() {
        rho[i] = rho0[i] + rho1[i];
    }
}

// ── Gate dispatch ─────────────────────────────────────────────────────────────

fn apply_op(rho: &mut [Complex64], n: usize, kind: GateKind, q: &[usize], p: &[f64]) -> Result<(), String> {
    match kind {
        GateKind::Id    => {}
        GateKind::X     => apply1_dm(rho, n, q[0], &gate_x()),
        GateKind::Y     => apply1_dm(rho, n, q[0], &gate_y()),
        GateKind::Z     => apply1_dm(rho, n, q[0], &gate_z()),
        GateKind::H     => apply1_dm(rho, n, q[0], &gate_h()),
        GateKind::S     => apply1_dm(rho, n, q[0], &gate_s()),
        GateKind::Sdg   => apply1_dm(rho, n, q[0], &gate_sdg()),
        GateKind::T     => apply1_dm(rho, n, q[0], &gate_t()),
        GateKind::Tdg   => apply1_dm(rho, n, q[0], &gate_tdg()),
        GateKind::Sx    => apply1_dm(rho, n, q[0], &gate_sx()),
        GateKind::Sxdg  => apply1_dm(rho, n, q[0], &gate_sxdg()),
        GateKind::Rx    => apply1_dm(rho, n, q[0], &gate_rx(p[0])),
        GateKind::Ry    => apply1_dm(rho, n, q[0], &gate_ry(p[0])),
        GateKind::Rz    => apply1_dm(rho, n, q[0], &gate_rz(p[0])),
        GateKind::Phase => apply1_dm(rho, n, q[0], &gate_phase(p[0])),
        GateKind::U     => apply1_dm(rho, n, q[0], &gate_u(p[0], p[1], p[2])),

        GateKind::Cx    => apply_ctrl1_dm(rho, n, q[0], q[1], &gate_x()),
        GateKind::Cy    => apply_ctrl1_dm(rho, n, q[0], q[1], &gate_y()),
        GateKind::Cz    => apply_ctrl1_dm(rho, n, q[0], q[1], &gate_z()),
        GateKind::Ch    => apply_ctrl1_dm(rho, n, q[0], q[1], &gate_h()),
        GateKind::Csx   => apply_ctrl1_dm(rho, n, q[0], q[1], &gate_sx()),
        GateKind::Crx   => apply_ctrl1_dm(rho, n, q[0], q[1], &gate_rx(p[0])),
        GateKind::Cry   => apply_ctrl1_dm(rho, n, q[0], q[1], &gate_ry(p[0])),
        GateKind::Crz   => apply_ctrl1_dm(rho, n, q[0], q[1], &gate_rz(p[0])),
        GateKind::Cp    => apply_ctrl1_dm(rho, n, q[0], q[1], &gate_phase(p[0])),
        GateKind::Cu    => apply_ctrl1_dm(rho, n, q[0], q[1], &gate_u(p[0], p[1], p[2])),
        GateKind::Swap  => apply_swap_dm(rho, n, q[0], q[1]),

        GateKind::Ccx   => apply_ccx_dm(rho, n, q[0], q[1], q[2]),
        GateKind::Cswap => apply_cswap_dm(rho, n, q[0], q[1], q[2]),
    }
    Ok(())
}

// ── Gate matrix definitions (inline copy, no pub dependency on statevector) ───

#[inline] fn c(re: f64, im: f64) -> Complex64 { Complex64::new(re, im) }

fn gate_x()   -> U2 { [[c(0.,0.),c(1.,0.)],[c(1.,0.),c(0.,0.)]] }
fn gate_y()   -> U2 { [[c(0.,0.),c(0.,-1.)],[c(0.,1.),c(0.,0.)]] }
fn gate_z()   -> U2 { [[c(1.,0.),c(0.,0.)],[c(0.,0.),c(-1.,0.)]] }
fn gate_h()   -> U2 { let f=FRAC_1_SQRT_2; [[c(f,0.),c(f,0.)],[c(f,0.),c(-f,0.)]] }
fn gate_s()   -> U2 { [[c(1.,0.),c(0.,0.)],[c(0.,0.),c(0.,1.)]] }
fn gate_sdg() -> U2 { [[c(1.,0.),c(0.,0.)],[c(0.,0.),c(0.,-1.)]] }
fn gate_t()   -> U2 { let e=Complex64::from_polar(1.,PI/4.); [[c(1.,0.),c(0.,0.)],[c(0.,0.),e]] }
fn gate_tdg() -> U2 { let e=Complex64::from_polar(1.,-PI/4.); [[c(1.,0.),c(0.,0.)],[c(0.,0.),e]] }
fn gate_sx()  -> U2 { [[c(0.5,0.5),c(0.5,-0.5)],[c(0.5,-0.5),c(0.5,0.5)]] }
fn gate_sxdg()-> U2 { [[c(0.5,-0.5),c(0.5,0.5)],[c(0.5,0.5),c(0.5,-0.5)]] }

fn gate_rx(t: f64) -> U2 {
    let cos=c((t/2.).cos(),0.); let isin=c(0.,-(t/2.).sin());
    [[cos,isin],[isin,cos]]
}
fn gate_ry(t: f64) -> U2 {
    let cos=c((t/2.).cos(),0.); let sin=c((t/2.).sin(),0.);
    [[cos,-sin],[sin,cos]]
}
fn gate_rz(t: f64) -> U2 {
    [[Complex64::from_polar(1.,-t/2.),c(0.,0.)],[c(0.,0.),Complex64::from_polar(1.,t/2.)]]
}
fn gate_phase(t: f64) -> U2 {
    [[c(1.,0.),c(0.,0.)],[c(0.,0.),Complex64::from_polar(1.,t)]]
}
fn gate_u(theta: f64, phi: f64, lambda: f64) -> U2 {
    let cos=(theta/2.).cos(); let sin=(theta/2.).sin();
    [[c(cos,0.),-Complex64::from_polar(sin,lambda)],
     [Complex64::from_polar(sin,phi),Complex64::from_polar(cos,phi+lambda)]]
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use cforge_core::Operation;
    use crate::{NativeStateVectorBackend, SimulationBackend};

    fn bell_circuit() -> Circuit {
        let mut c = Circuit::new(2);
        c.push(Operation::new(GateKind::H,  vec![0], vec![]));
        c.push(Operation::new(GateKind::Cx, vec![0,1], vec![]));
        c
    }

    #[test]
    fn zero_noise_bell_probabilities_correct() {
        let backend = DensityMatrixBackend::default();
        let r = backend.run(&bell_circuit(), 0, 0).unwrap();
        let probs: Vec<f64> = r.statevector.iter().map(|a| a.norm_sqr()).collect();
        // Bell state: |00⟩ and |11⟩ each ~0.5, others ~0
        assert!((probs[0] - 0.5).abs() < 1e-9, "p(00)={}", probs[0]);
        assert!((probs[3] - 0.5).abs() < 1e-9, "p(11)={}", probs[3]);
        assert!(probs[1] < 1e-9 && probs[2] < 1e-9);
    }

    #[test]
    fn depolarizing_noise_degrades_fidelity() {
        use crate::noise::NoisyConfig;
        let noisy = DensityMatrixBackend { config: Some(NoisyConfig::depolarizing(0.3, 0.3)) };
        let r = noisy.run(&bell_circuit(), 0, 0).unwrap();
        let probs: Vec<f64> = r.statevector.iter().map(|a| a.norm_sqr()).collect();
        // With heavy noise, distribution should be less extreme than 0.5/0.5
        let max_p = probs.iter().cloned().fold(0.0_f64, f64::max);
        assert!(max_p < 0.5, "noisy Bell state max probability should be < 0.5, got {max_p}");
    }

    #[test]
    fn amplitude_damping_decays_excited_state() {
        use crate::noise::NoisyConfig;
        let mut c = Circuit::new(1);
        c.push(Operation::new(GateKind::X, vec![0], vec![])); // prepare |1⟩
        let noisy = DensityMatrixBackend { config: Some(NoisyConfig::amplitude_damping(0.9)) };
        let r = noisy.run(&c, 0, 0).unwrap();
        let p1 = r.statevector[1].norm_sqr(); // prob of |1⟩
        // After γ=0.9 damping, most population should have decayed to |0⟩
        assert!(p1 < 0.2, "p(|1⟩) after strong damping = {p1}");
    }

    #[test]
    fn density_matrix_matches_statevector_noiseless() {
        // For a pure state |ψ⟩, ρ = |ψ⟩⟨ψ| → diagonal = |ψ_i|².
        // DM backend returns sqrt(diag) so |amplitude|² == native |amplitude|².
        let mut c = Circuit::new(2);
        c.push(Operation::new(GateKind::H,  vec![0], vec![]));
        c.push(Operation::new(GateKind::Rz, vec![0], vec![std::f64::consts::FRAC_PI_3]));
        c.push(Operation::new(GateKind::Cx, vec![0,1], vec![]));

        let native = NativeStateVectorBackend.run(&c, 0, 0).unwrap();
        let dm     = DensityMatrixBackend::default().run(&c, 0, 0).unwrap();

        for (i, (n, d)) in native.statevector.iter().zip(&dm.statevector).enumerate() {
            let p_native = n.norm_sqr();
            let p_dm     = d.norm_sqr();
            assert!((p_native - p_dm).abs() < 1e-9,
                    "index {i}: native={p_native:.6} dm={p_dm:.6}");
        }
    }

    #[test]
    fn shots_deterministic() {
        let backend = DensityMatrixBackend::default();
        let c = bell_circuit();
        let r1 = backend.run(&c, 1000, 42).unwrap();
        let r2 = backend.run(&c, 1000, 42).unwrap();
        assert_eq!(r1.counts, r2.counts);
    }

    #[test]
    fn density_matrix_backend_name() {
        assert_eq!(DensityMatrixBackend::default().name(), "density-matrix");
    }
}
