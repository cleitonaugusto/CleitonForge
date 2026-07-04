// Quantum Fourier Transform on 3 qubits, input state |100>
// Expected: uniform superposition (all 8 states with prob 1/8)
OPENQASM 2.0;
qreg q[3];
x q[0];
h q[0];
cp(pi/2) q[1],q[0];
cp(pi/4) q[2],q[0];
h q[1];
cp(pi/2) q[2],q[1];
h q[2];
swap q[0],q[2];
