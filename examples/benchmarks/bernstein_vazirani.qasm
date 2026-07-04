// Bernstein-Vazirani: hidden string s=101
// Phase oracle: Z(q0) Z(q2) gives phase (-1)^{x0+x2} = (-1)^{s.x}
// Expected result: |101> with 100% probability
OPENQASM 2.0;
qreg q[3];
h q[0];
h q[1];
h q[2];
z q[0];
z q[2];
h q[0];
h q[1];
h q[2];
