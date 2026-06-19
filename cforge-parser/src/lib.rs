//! Multi-format quantum circuit parsing for CleitonForge.
//!
//! Entry points:
//! - [`qasm2::parse_qasm2`] — OpenQASM 2.0 source string → [`Circuit`]
//! - [`qasm3::parse_qasm3`] — OpenQASM 3.0 source string → [`Circuit`]
//!
//! All parsers produce CleitonForge's canonical [`Circuit`] IR so that
//! backends and metrics never need to understand a source format.

pub mod error;
pub mod qasm2;
pub mod qasm3;
mod translate;

pub use cforge_core::{Circuit, GateKind, Operation};
pub use error::ParseError;
pub use qasm2::parse_qasm2;
pub use qasm3::parse_qasm3;
