OPENQASM 2.0;

// QAOA p=1 — MaxCut on 4-node ring graph
// Graph edges: (0,1), (1,2), (2,3), (3,0)
// Parameters: gamma = pi/4 (phase), beta = pi/8 (mixer)
// These values are analytically optimal for p=1 on a 4-ring.
qreg q[4];

// Initial state |+>^4
h q[0];
h q[1];
h q[2];
h q[3];

// Phase separator: e^{-i*gamma*Z_i Z_j} for each edge
// Implemented as: CX(i,j) RZ(2*gamma, j) CX(i,j)
// gamma = pi/4  =>  2*gamma = pi/2 = 1.5707963267948966

// Edge (0,1)
cx q[0],q[1];
rz(1.5707963267948966) q[1];
cx q[0],q[1];

// Edge (1,2)
cx q[1],q[2];
rz(1.5707963267948966) q[2];
cx q[1],q[2];

// Edge (2,3)
cx q[2],q[3];
rz(1.5707963267948966) q[3];
cx q[2],q[3];

// Edge (3,0)
cx q[3],q[0];
rz(1.5707963267948966) q[0];
cx q[3],q[0];

// Mixer: e^{-i*beta*X_i} = RX(2*beta) for each qubit
// beta = pi/8  =>  2*beta = pi/4 = 0.7853981633974483
rx(0.7853981633974483) q[0];
rx(0.7853981633974483) q[1];
rx(0.7853981633974483) q[2];
rx(0.7853981633974483) q[3];
