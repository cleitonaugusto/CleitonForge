# research/

Experiments behind the claim that a differential-testing oracle must declare its
own power. Each script is self-contained and prints its numbers; nothing here is
typed in by hand.

| script | what it answers |
|---|---|
| `oracle_blindness_theory.py` | For which gates does inverting the Rz sign conjugate the circuit? Gives the exact blind family. |
| `oracle_blindness_validate.py` | Does the predicted blind family actually show zero counts divergence over random circuits? |
| `oracle_power_hierarchy.py` | Power of the counts, state and operator oracles, with a Wilson interval. |

## Result, 2026-09-01

Fault injected: invert the sign of `Rz`. Nothing else changes.

**Blind family** — counts oracle has provably zero power over
`{I, H, X, Z, Ry, Rz, CX, CZ, SWAP}`. Adding `Y` keeps it blind, because
`conj(Y) = -Y` is a global phase while `Y` is uncontrolled. Adding `S`, `T`,
`SX`, `Rx` or `P` breaks blindness.

**Power**, 600 random 4-qubit circuits of depth 14 per family:

| family | counts | state | operator |
|---|---|---|---|
| `{H,X,Z,Ry,Rz,CX,CZ}` | **0.0%** | 50.8% | 89.0% |
| `+ T` | 2.3% | 44.2% | 81.3% |
| `+ Rx` | 12.8% | 58.5% | 86.3% |

Wilson 95% interval on the zero: **[0.00%, 0.64%]**.

Three things follow, and the second and third were not expected.

1. The analytic prediction holds exactly. Zero, not nearly zero.
2. **Leaving the blind family does not buy power.** Adding `T` moves the counts
   oracle from 0% to 2.3%; adding `Rx`, to 12.8%. Breaking blindness and gaining
   detection are different things, so the problem is not a pathological family —
   it is the oracle.
3. **No oracle here is complete.** The operator oracle sits at 81–89%, missing
   the cases where the inversion yields an operator equal up to global phase.

Point 3 strengthens the thesis rather than weakening it: if no oracle has full
power, declaring the power you do have stops being good practice and becomes
necessary to interpret any negative result.
