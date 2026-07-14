OPENQASM 2.0;

// Bernstein-Vazirani — 4 input qubits + 1 ancilla
// Secret string s = 1011 (binary) = bits 0,1,3 set
// Algorithm recovers s in a single oracle query.
qreg q[4];
qreg anc[1];

// Ancilla in |-> = H|1>
x anc[0];
h anc[0];

// Input qubits in |+>^4
h q[0];
h q[1];
h q[2];
h q[3];

// Oracle: for each bit i where s_i=1, apply CX q[i] -> anc
cx q[0],anc[0];
cx q[1],anc[0];
// q[2] skipped — s_2 = 0
cx q[3],anc[0];

// Hadamard to extract s
h q[0];
h q[1];
h q[2];
h q[3];
// Expected statevector peak: |1011> (q3 q2 q1 q0 in q0=LSB order)
