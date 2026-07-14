OPENQASM 2.0;

// Quantum Fourier Transform — 4 qubits
// Convention: QFT|x> = (1/4) sum_k e^{2pi i x k / 16} |k>
// Bit-reversal swaps at the end produce the standard ordered output.
qreg q[4];

// QFT on q[0]
h q[0];
cp(1.5707963267948966) q[1],q[0];
cp(0.7853981633974483) q[2],q[0];
cp(0.3926990816987242) q[3],q[0];

// QFT on q[1]
h q[1];
cp(1.5707963267948966) q[2],q[1];
cp(0.7853981633974483) q[3],q[1];

// QFT on q[2]
h q[2];
cp(1.5707963267948966) q[3],q[2];

// QFT on q[3]
h q[3];

// Bit-reversal
swap q[0],q[3];
swap q[1],q[2];
