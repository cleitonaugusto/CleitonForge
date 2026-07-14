OPENQASM 2.0;

// Grover's algorithm — 2 qubits, oracle marks |11>
// One Grover iteration is optimal for n=2 (amplitude ~1.0 on |11>)
qreg q[2];

// Step 1: uniform superposition
h q[0];
h q[1];

// Step 2: phase oracle for |11> (CZ = controlled-Z)
cz q[0],q[1];

// Step 3: diffusion operator  H X CZ X H  on both qubits
h q[0];
h q[1];
x q[0];
x q[1];
cz q[0],q[1];
x q[0];
x q[1];
h q[0];
h q[1];
