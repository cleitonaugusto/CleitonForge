"""Shared primitives for the oracle-power experiments.

A circuit is a list of ops, not a matrix. That matters: mutations act on the
circuit, the way a compiler bug does, instead of being a flag threaded through
gate construction. Simulation is a separate step, so the same circuit can be
run through several oracles without rebuilding it.

Single-qubit gates are applied to the full DIM x DIM operator by tensor
contraction, O(4^n) per gate rather than the O(8^n) of a full matrix product.
"""
import numpy as np

# ---------------------------------------------------------------- gates
I2 = np.eye(2, dtype=complex)
H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]])
Z = np.diag([1, -1]).astype(complex)
S = np.diag([1, 1j])
SDG = np.diag([1, -1j])
T = np.diag([1, np.exp(1j * np.pi / 4)])
TDG = np.diag([1, np.exp(-1j * np.pi / 4)])
SX = 0.5 * np.array([[1 + 1j, 1 - 1j], [1 - 1j, 1 + 1j]])


def rx(t):
    c, s = np.cos(t / 2), np.sin(t / 2)
    return np.array([[c, -1j * s], [-1j * s, c]])


def ry(t):
    c, s = np.cos(t / 2), np.sin(t / 2)
    return np.array([[c, -s], [s, c]], dtype=complex)


def rz(t):
    return np.diag([np.exp(-1j * t / 2), np.exp(1j * t / 2)])


ONE_Q = {"H": H, "X": X, "Y": Y, "Z": Z, "S": S, "T": T, "SX": SX}
ROT = {"Rx": rx, "Ry": ry, "Rz": rz}


# ---------------------------------------------------------------- circuits
# An op is a tuple:
#   ("1q",   name, qubit, angle_or_None)
#   ("2q",   name, control, target)
def random_circuit(rng, n, depth, gates):
    """Build a random circuit over `gates` on `n` qubits."""
    ops = []
    for _ in range(depth):
        g = gates[rng.integers(len(gates))]
        if g in ("CX", "CZ"):
            if n < 2:
                continue
            a, b = rng.choice(n, 2, replace=False)
            ops.append(("2q", g, int(a), int(b)))
        elif g in ROT:
            ops.append(("1q", g, int(rng.integers(n)),
                        float(rng.uniform(-np.pi, np.pi))))
        else:
            ops.append(("1q", g, int(rng.integers(n)), None))
    return ops


def _apply_1q(U, m, q, n):
    """U <- (I x .. x m x .. x I) @ U, by contraction on axis q."""
    dim = U.shape[0]
    U = U.reshape((2,) * n + (dim,))
    U = np.moveaxis(U, q, 0).reshape(2, -1)
    U = m @ U
    U = U.reshape((2,) + (2,) * (n - 1) + (dim,))
    return np.moveaxis(U, 0, q).reshape(dim, dim)


def _apply_2q(U, name, c, t, n):
    dim = U.shape[0]
    idx = np.arange(dim)
    cbit = (idx >> (n - 1 - c)) & 1
    if name == "CX":
        tgt = idx ^ (cbit << (n - 1 - t))
        return U[tgt, :]
    phase = np.where(cbit & ((idx >> (n - 1 - t)) & 1), -1.0, 1.0)
    return U * phase[:, None]


def simulate(ops, n):
    """Full unitary of the circuit, as a DIM x DIM array."""
    dim = 2 ** n
    U = np.eye(dim, dtype=complex)
    for op in ops:
        if op[0] == "1q":
            _, name, q, ang = op
            m = ROT[name](ang) if name in ROT else ONE_Q[name]
            U = _apply_1q(U, m, q, n)
        elif op[0] == "1qm":          # gate carrying its own matrix
            _, m, q, _ = op
            U = _apply_1q(U, m, q, n)
        else:
            _, name, c, t = op
            U = _apply_2q(U, name, c, t, n)
    return U


# ---------------------------------------------------------------- mutations
# Each returns a NEW op list, or None when it does not apply to this circuit.
def mut_drop_gate(ops, rng):
    if not ops:
        return None
    i = int(rng.integers(len(ops)))
    return ops[:i] + ops[i + 1:]


def mut_flip_rotation_sign(ops, rng):
    cand = [i for i, o in enumerate(ops) if o[0] == "1q" and o[1] in ROT]
    if not cand:
        return None
    i = int(rng.choice(cand))
    _, name, q, ang = ops[i]
    return ops[:i] + [("1q", name, q, -ang)] + ops[i + 1:]


def mut_flip_all_rz_signs(ops, rng):
    """The convention bug: every Rz built with the opposite sign.

    Distinct from flipping one gate, and the distinction is the whole point.
    Inverting every Rz conjugates a circuit made of otherwise-real gates, and
    |z|^2 = |conj(z)|^2, so no measurement probability moves. Inverting a
    single Rz conjugates nothing and is ordinarily visible. Calling both
    "rotation sign" hides a theorem.
    """
    cand = [i for i, o in enumerate(ops) if o[0] == "1q" and o[1] == "Rz"]
    if not cand:
        return None
    out = list(ops)
    for i in cand:
        _, name, q, ang = out[i]
        out[i] = ("1q", name, q, -ang)
    return out


def mut_swap_qubit_order(ops, rng):
    cand = [i for i, o in enumerate(ops) if o[0] == "2q"]
    if not cand:
        return None
    i = int(rng.choice(cand))
    _, name, c, t = ops[i]
    return ops[:i] + [("2q", name, t, c)] + ops[i + 1:]


def mut_adjoint_neighbour(ops, rng):
    """S <-> Sdg, T <-> Tdg: the classic off-by-a-dagger."""
    swap = {"S": SDG, "T": TDG}
    cand = [i for i, o in enumerate(ops) if o[0] == "1q" and o[1] in swap]
    if not cand:
        return None
    i = int(rng.choice(cand))
    _, name, q, _ = ops[i]
    return ops[:i] + [("1qm", swap[name], q, None)] + ops[i + 1:]


def mut_truncate_angle(ops, rng):
    """Round an angle to the nearest multiple of pi/2, as a lossy pass would."""
    cand = [i for i, o in enumerate(ops) if o[0] == "1q" and o[1] in ROT]
    if not cand:
        return None
    i = int(rng.choice(cand))
    _, name, q, ang = ops[i]
    snapped = float(np.round(ang / (np.pi / 2)) * (np.pi / 2))
    if abs(snapped - ang) < 1e-12:
        return None
    return ops[:i] + [("1q", name, q, snapped)] + ops[i + 1:]


MUTATIONS = {
    "drop gate": mut_drop_gate,
    "rotation sign (one gate)": mut_flip_rotation_sign,
    "Rz convention (all gates)": mut_flip_all_rz_signs,
    "qubit order": mut_swap_qubit_order,
    "adjoint neighbour": mut_adjoint_neighbour,
    "angle truncation": mut_truncate_angle,
}


# ---------------------------------------------------------------- oracles
def probs(U):
    """Outcome probabilities from |0...0>."""
    return np.abs(U[:, 0]) ** 2


def oracle_sampled_counts(Ua, Ub, rng, shots=1024, alpha=0.01):
    """Weakest oracle, and the realistic one: finite shots, chi-square style.

    Included as a negative control. Real campaigns compare histograms from a
    few thousand shots; this shows what that costs relative to exact counts.
    """
    pa, pb = probs(Ua), probs(Ub)
    sa = rng.multinomial(shots, pa / pa.sum())
    sb = rng.multinomial(shots, pb / pb.sum())
    tot = sa + sb
    keep = tot > 0
    exp = tot[keep] / 2.0
    chi2 = (((sa[keep] - exp) ** 2 + (sb[keep] - exp) ** 2) / exp).sum()
    dof = max(int(keep.sum()) - 1, 1)
    # 99% cut of chi2(dof), Wilson-Hilferty
    z = 2.3263
    cut = dof * (1 - 2 / (9 * dof) + z * np.sqrt(2 / (9 * dof))) ** 3
    return chi2 > cut


def oracle_exact_counts(Ua, Ub, tol=1e-10):
    return not np.allclose(probs(Ua), probs(Ub), atol=tol)


def oracle_state(Ua, Ub, tol=1e-10):
    a, b = Ua[:, 0], Ub[:, 0]
    return abs(abs(np.vdot(a, b)) ** 2 - 1) > tol


def oracle_operator(Ua, Ub, tol=1e-9):
    i = np.unravel_index(np.argmax(np.abs(Ub)), Ub.shape)
    if abs(Ub[i]) < tol:
        return not np.allclose(Ua, Ub, atol=tol)
    return not np.allclose(Ua, (Ua[i] / Ub[i]) * Ub, atol=tol)


# ---------------------------------------------------------------- statistics
def wilson(k, n, z=1.96):
    """Wilson score interval. Correct at k=0, unlike the normal approximation."""
    if n == 0:
        return (0.0, 0.0, 0.0)
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (p, max(0.0, centre - half), min(1.0, centre + half))
