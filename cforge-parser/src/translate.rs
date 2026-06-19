//! Mapping from gate name strings to cforge-core `GateKind`.

use cforge_core::GateKind;

/// Looks up a gate name (case-insensitive) and returns the matching `GateKind`.
/// Returns `None` for gates not in the stdgates.inc universal set.
pub fn gate_kind_from_name(name: &str) -> Option<GateKind> {
    match name.to_lowercase().as_str() {
        "id" | "i" => Some(GateKind::Id),
        "x" => Some(GateKind::X),
        "y" => Some(GateKind::Y),
        "z" => Some(GateKind::Z),
        "h" => Some(GateKind::H),
        "s" => Some(GateKind::S),
        "sdg" => Some(GateKind::Sdg),
        "t" => Some(GateKind::T),
        "tdg" => Some(GateKind::Tdg),
        "sx" => Some(GateKind::Sx),
        "sxdg" => Some(GateKind::Sxdg),
        "rx" => Some(GateKind::Rx),
        "ry" => Some(GateKind::Ry),
        "rz" => Some(GateKind::Rz),
        "p" | "phase" | "u1" => Some(GateKind::Phase),
        "u" | "u3" => Some(GateKind::U),
        "u2" => None, // u2(phi,lam) = u(pi/2,phi,lam); handled by callers if needed
        "cx" | "cnot" => Some(GateKind::Cx),
        "cy" => Some(GateKind::Cy),
        "cz" => Some(GateKind::Cz),
        "ch" => Some(GateKind::Ch),
        "csx" => Some(GateKind::Csx),
        "crx" => Some(GateKind::Crx),
        "cry" => Some(GateKind::Cry),
        "crz" => Some(GateKind::Crz),
        "cp" | "cphase" => Some(GateKind::Cp),
        "cu" | "cu3" => Some(GateKind::Cu),
        "swap" => Some(GateKind::Swap),
        "ccx" | "toffoli" | "ccnot" => Some(GateKind::Ccx),
        "cswap" | "fredkin" => Some(GateKind::Cswap),
        _ => None,
    }
}
