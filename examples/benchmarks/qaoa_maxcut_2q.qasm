// QAOA p=1 for MaxCut on a single edge (q0-q1)
// Optimal angles: gamma=-3pi/4, beta=-pi/8 (statevector convention)
// Expected: 100% cut states |01> + |10> on statevector backend
// Note: quantrs2 backend gives different result due to Rz sign convention
//       (documented inter-framework divergence, not a bug)
OPENQASM 2.0;
qreg q[2];
h q[0]; h q[1];
cx q[0],q[1];
rz(-3*pi/2) q[1];
cx q[0],q[1];
rx(-pi/4) q[0];
rx(-pi/4) q[1];
