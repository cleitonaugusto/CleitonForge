#!/usr/bin/env python3
"""The exact-operator cross-check has to reject QCEC's misfire and nothing else.

confirm_bug exists because QCEC reports not_equivalent when a small-angle
rotation is compared against the same operator packed into a UnitaryGate, which
is what ConsolidateBlocks emits. That misfire has a true error of exactly 0.0.

The trap is using a loose tolerance to absorb it. Operator.equiv defaults to
rtol=1e-5, which also absorbs genuine divergences of that size, and a genuine
divergence discarded as an oracle false positive never reaches a human at all.
bug-zoo/qiskit-cp-small-angle-002 is that exact shape.

Run: python3 tools/test_confirm_bug.py
"""
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])

PI_2 = 1.5707963267948966


def check(label: str, condition: bool) -> bool:
    print(f"  [{'OK' if condition else 'FALHA'}] {label}")
    return condition


def main() -> int:
    try:
        from qiskit import QuantumCircuit
        from qiskit.quantum_info import Operator
    except ImportError:
        print("qiskit ausente, pulando")
        return 0

    from oracle_qcec import confirm_bug

    ok = True

    print("rejeita a falha do QCEC (erro verdadeiro exatamente zero)")
    for theta in (1e-4, 5e-5, 1e-5):
        a = QuantumCircuit(2)
        a.crz(theta, 0, 1)
        b = QuantumCircuit(2)
        b.unitary(Operator(a).data, [0, 1])
        ok = check(f"crz({theta:g}) vs a propria matriz", confirm_bug(a, b) is False) and ok

    print("mas nao engole divergencia real da mesma ordem de grandeza")
    for delta in (1e-3, 1e-4, 1e-5, 1e-6):
        c = QuantumCircuit(1)
        c.rz(PI_2, 0)
        d = QuantumCircuit(1)
        d.rz(PI_2 + delta, 0)
        ok = check(f"rz(pi/2) vs rz(pi/2 + {delta:g})", confirm_bug(c, d) is True) and ok

    print("circuito largo demais devolve None, nao False")
    wide = QuantumCircuit(30)
    wide.h(0)
    ok = check("None quando nao da para checar", confirm_bug(wide, wide) is None) and ok

    print("OK" if ok else "FALHOU")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
