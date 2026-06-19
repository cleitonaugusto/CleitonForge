/// An error that occurred while parsing a quantum circuit source file.
#[derive(Debug, Clone, PartialEq)]
pub enum ParseError {
    /// The source text is syntactically invalid.
    SyntaxError(String),
    /// A qubit name referenced in a gate application is not declared.
    UndeclaredQubit(String),
    /// A gate name could not be mapped to a known `GateKind`.
    UnknownGate(String),
    /// A numeric parameter expression could not be evaluated.
    InvalidParam(String),
}

impl std::fmt::Display for ParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParseError::SyntaxError(msg) => write!(f, "syntax error: {msg}"),
            ParseError::UndeclaredQubit(name) => write!(f, "undeclared qubit register: {name}"),
            ParseError::UnknownGate(name) => write!(f, "unknown gate: {name}"),
            ParseError::InvalidParam(expr) => write!(f, "invalid parameter expression: {expr}"),
        }
    }
}

impl std::error::Error for ParseError {}
