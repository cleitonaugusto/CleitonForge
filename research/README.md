# research/

Experiments behind the claim that a differential-testing oracle must declare
its own power. Every number below is printed by a script in this directory.

| script | what it answers |
|---|---|
| `qcore.py` | Shared primitives: circuits as op lists, simulation, mutations, oracles, Wilson intervals. Not a script. |
| `oracle_blindness_theory.py` | For which gates does inverting the Rz sign conjugate the circuit? Gives the exact blind family. |
| `oracle_blindness_validate.py` | Does the predicted blind family actually show zero counts divergence over random circuits? |
| `oracle_power_hierarchy.py` | The hierarchy counts < state < operator for the Rz convention fault, in the smallest setting where blindness is true. |
| `oracle_power_study.py` | Factorial study: 4 families x 3 widths x 3 depths x 6 fault classes x 4 oracles. |

## Method

**Power is measured against non-equivalent mutants.** A mutation that leaves the
operator unchanged up to global phase cannot be detected by any correct oracle,
so counting it as a miss understates every oracle at once. Equivalence is decided
by the operator oracle, which is therefore the *reference*: its power is 1 by
construction, not a measured result. What the study measures is how far the
cheaper oracles fall below it.

**The sampled oracle is scored for false positives.** Equivalent mutants are a
free negative control: the operators agree, so any alarm is wrong. Measured at
231/26666 = **0.9%** [95% CI 0.8%, 1.0%], matching the nominal 1% of the
chi-square cut. A raw detection rate that includes false positives is not a
power.

## Results, 2026-09-01

284,595 mutants over 216 cells; 26,666 equivalent (9.4%); 257,929 scored.

### 1. The blindness theorem, at scale

Fault: **every** `Rz` built with the opposite sign — a convention bug.

| family | counts | state |
|---|---|---|
| `{H,X,Z,Ry,Rz,CX,CZ}` | **0/16418 = 0.0%**, CI [0.000%, 0.023%] | 72.8% |
| `{Rx,Ry,Rz,CX}` (has `Rx`) | 65.4% | 86.0% |

Zero, not nearly zero, across three widths and three depths. The state oracle
flags the same mutants, which is what proves they are not equivalent.

### 2. Inverting every Rz is not the same as inverting one

| fault | counts | state |
|---|---|---|
| `Rz` convention (all gates) | 33.9% | 82.8% |
| rotation sign (one gate) | 52.0% | 81.1% |

Blindness is a property of the *global* error. Inverting every `Rz` conjugates a
circuit made of otherwise real gates, and `|z|^2 = |conj(z)|^2`, so no
probability moves. Inverting a single `Rz` conjugates nothing. Calling both
"rotation sign" hides the theorem.

### 3. A second blindness, structural, in Clifford circuits

Fault: one `S` replaced by `S-dagger`. Counts power, exact:

| family | n | d=10 | d=20 | d=40 |
|---|---|---|---|---|
| Clifford | 3 | 0.9% | 6.1% | 14.6% |
| Clifford | 5 | 0.3% | 1.9% | 8.8% |
| Clifford+T | 5 | 0.3% | 2.1% | 10.1% |

Near-total blindness, and **not** the conjugation theorem: that hypothesis
predicts 80-100% power here, two orders of magnitude off. Measured directly
instead:

- 100% of these circuits produce a flat distribution (a stabilizer state gives
  uniform probability over its support)
- 96% keep **identical** counts after the mutation

A phase-only fault does not move the support of a stabilizer state, so counts
see nothing. This is about the circuit family, not about a gate, and it applies
to any campaign whose targets are Clifford — which includes error-correction
tooling, where detector samples are the usual oracle.

### 4. Power moves with depth and width, monotonically

| | counts | state |
|---|---|---|
| depth 10 | 34.2% | 54.5% |
| depth 20 | 46.1% | 69.1% |
| depth 40 | 58.8% | 81.8% |

| | counts | state |
|---|---|---|
| 3 qubits | 50.4% | 72.7% |
| 4 qubits | 46.3% | 68.4% |
| 5 qubits | 42.8% | 64.8% |

An oracle's power is not a property of the oracle. It is a property of the pair
(oracle, workload), which is why it has to be reported per campaign rather than
once.

### 5. Fault classes are not interchangeable

| fault | equivalent | counts | state |
|---|---|---|---|
| qubit order | 37.4% | 63.0% | 69.4% |
| drop gate | 0.0% | 47.7% | 62.3% |
| angle truncation | 0.0% | 69.8% | 81.2% |
| adjoint neighbour | 0.0% | 6.3% | 41.7% |

`qubit order` produces 37.4% equivalent mutants because `CZ` is symmetric in its
two qubits — swapping them is a no-op. That is a sanity check the harness must
pass, and a reminder that "detected / generated" is the wrong ratio.

## Correction to the previous version of this file

The earlier README claimed:

> "No oracle here is complete. The operator oracle sits at 81-89%, missing the
> cases where the inversion yields an operator equal up to global phase."

**That was wrong.** The misses were circuits containing no `Rz` at all, where
inverting the `Rz` sign changes nothing and the operator is correctly identical.
With depth 14 over 7 gates, `(6/7)^14 = 11.6%`, predicting 88.4% against the
89.0% published. Within noise at n=600.

The cause was in the design, not the arithmetic: the fault was a flag threaded
through gate construction (`rz(t, bad=True)`) rather than a mutation applied to
a circuit, so there was no way to ask whether the mutation changed anything. A
circuit with no `Rz` and an oracle that failed counted the same.

Separating circuit from simulation made the question appear on its own.
