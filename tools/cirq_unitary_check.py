"""Compara OPERADOR, nao estado. Complementa o --check do harness do Cirq.

Comparar o estado a partir de |0...0> so exercita uma coluna do unitario. Medido
por mutacao sobre 2000 circuitos: remover um gate e detectado em 65% dos casos, e
os 35% restantes sao gates que caem como identidade no estado alcancado ate ali.
Um defeito que so aparece a partir de outra entrada passa despercebido.

Roda so nos casos em que o export manteve todos os qubits (96% da amostra), para
nao precisar embutir a matriz no espaco completo.

Uso: python3 unit_check.py casos.jsonl
"""
import json
import re
import sys

import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import Operator


def qubit_map(qasm):
    m = re.search(r"^// Qubits: \[(.*)\]", qasm, re.M)
    return [int(x) for x in re.findall(r"q\((\d+)\)", m.group(1))] if m else None


def as_matrix(pairs):
    """Reconstroi a matriz complexa gravada como pares [re, im]."""
    return np.array([[complex(re_, im_) for re_, im_ in row] for row in pairs])


def diff_op(a, b):
    """Maior diferenca entre entradas, a menos de fase global."""
    tr = np.vdot(a.reshape(-1), b.reshape(-1))
    ph = tr / abs(tr) if abs(tr) > 1e-12 else 1.0
    return float(np.max(np.abs(a - b * np.conj(ph))))


def rebuild_operator(qasm, num_qubits):
    qc = QuantumCircuit.from_qasm_str(qasm)
    mp = qubit_map(qasm)
    if mp is None or len(mp) != num_qubits or qc.num_qubits != num_qubits:
        return None
    return Operator(qc).reverse_qargs().data


def main(path):
    checked = skipped = bad = 0
    worst = 0.0
    for line in open(path):
        r = json.loads(line)
        if r.get("calibration") or "error" in r or "unitary" not in r:
            skipped += 1
            continue
        got = rebuild_operator(r["qasm2"], r["num_qubits"])
        if got is None:
            skipped += 1
            continue
        want = as_matrix(r["unitary"])
        d = diff_op(want, got)
        worst = max(worst, d)
        checked += 1
        if d > 1e-6:
            bad += 1
            print(f"  [MISMATCH] caso {r['i']} (seed {r['seed']}), diferenca {d:.3e}")
            for n_, q_, p_ in r["ops"]:
                print(f"      {n_}({', '.join(f'{x:.6f}' for x in p_)}) em {q_}")
    print(f"\noperadores comparados: {checked}, divergencias: {bad}, pulados: {skipped}")
    print(f"maior diferenca observada: {worst:.3e}")
    return 1 if (bad or not checked) else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1]))
