// Grover's algorithm: 3 qubits, target state |101>, 2 iterations
// Oracle: CCZ via phase kickback (H-CCX-H). Expected: ~94.5% |101>
// Reference: sin^2((2*2+1)*arcsin(1/sqrt(8))) = 0.9453
OPENQASM 2.0;
qreg q[3];
h q[0]; h q[1]; h q[2];
// iteration 1
x q[1]; h q[2]; ccx q[0],q[1],q[2]; h q[2]; x q[1];
h q[0]; h q[1]; h q[2];
x q[0]; x q[1]; x q[2];
h q[2]; ccx q[0],q[1],q[2]; h q[2];
x q[0]; x q[1]; x q[2];
h q[0]; h q[1]; h q[2];
// iteration 2
x q[1]; h q[2]; ccx q[0],q[1],q[2]; h q[2]; x q[1];
h q[0]; h q[1]; h q[2];
x q[0]; x q[1]; x q[2];
h q[2]; ccx q[0],q[1],q[2]; h q[2];
x q[0]; x q[1]; x q[2];
h q[0]; h q[1]; h q[2];
