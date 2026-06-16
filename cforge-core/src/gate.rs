//! Universal gate set, matching the minimal gate set defined by
//! OpenQASM 3's `stdgates.inc`.

/// The kind of a quantum gate, independent of which qubits it acts on or
/// which framework it originated from.
///
/// Each variant documents how many qubits it acts on and how many
/// numeric parameters it expects, via [`GateKind::num_qubits`] and
/// [`GateKind::num_params`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum GateKind {
    Id,
    X,
    Y,
    Z,
    H,
    S,
    Sdg,
    T,
    Tdg,
    Sx,
    Sxdg,
    Rx,
    Ry,
    Rz,
    Phase,
    /// General single-qubit unitary `U(theta, phi, lambda)`.
    U,
    Cx,
    Cy,
    Cz,
    Ch,
    Csx,
    Crx,
    Cry,
    Crz,
    Cp,
    /// Controlled general single-qubit unitary `CU(theta, phi, lambda, gamma)`.
    Cu,
    Swap,
    Ccx,
    Cswap,
}

impl GateKind {
    /// Number of qubits this gate kind acts on.
    pub fn num_qubits(&self) -> usize {
        use GateKind::*;
        match self {
            Id | X | Y | Z | H | S | Sdg | T | Tdg | Sx | Sxdg | Rx | Ry | Rz | Phase | U => 1,
            Cx | Cy | Cz | Ch | Csx | Crx | Cry | Crz | Cp | Cu | Swap => 2,
            Ccx | Cswap => 3,
        }
    }

    /// Number of numeric parameters this gate kind expects.
    pub fn num_params(&self) -> usize {
        use GateKind::*;
        match self {
            Rx | Ry | Rz | Phase | Crx | Cry | Crz | Cp => 1,
            U => 3,
            Cu => 4,
            _ => 0,
        }
    }

    /// Canonical OpenQASM 3 `stdgates.inc` name for this gate kind.
    pub fn qasm_name(&self) -> &'static str {
        use GateKind::*;
        match self {
            Id => "id",
            X => "x",
            Y => "y",
            Z => "z",
            H => "h",
            S => "s",
            Sdg => "sdg",
            T => "t",
            Tdg => "tdg",
            Sx => "sx",
            Sxdg => "sxdg",
            Rx => "rx",
            Ry => "ry",
            Rz => "rz",
            Phase => "p",
            U => "u",
            Cx => "cx",
            Cy => "cy",
            Cz => "cz",
            Ch => "ch",
            Csx => "csx",
            Crx => "crx",
            Cry => "cry",
            Crz => "crz",
            Cp => "cp",
            Cu => "cu",
            Swap => "swap",
            Ccx => "ccx",
            Cswap => "cswap",
        }
    }
}
