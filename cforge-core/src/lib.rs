//! Canonical intermediate representation and shared types for CleitonForge.
//!
//! Every parser translates into [`ir::Circuit`], and every backend and
//! metric operates only on that IR — never on a framework-native
//! representation. That's what makes results comparable across
//! frameworks.

pub mod gate;
pub mod ir;
pub mod metrics;

pub use gate::GateKind;
pub use ir::{Circuit, CircuitError, Operation, Qubit};
pub use metrics::MetricsResult;
