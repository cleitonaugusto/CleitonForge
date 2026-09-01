"""Qual e exatamente a familia cega para o oraculo de contagens?

Hipotese analitica: inverter o sinal do Rz conjuga o unitario inteiro se, e
somente se, TODO gate do circuito satisfaz uma de duas condicoes:
   conj(g) = g            (o gate e real)
   conj(g) = g_invertido  (o backend defeituoso ja produz o conjugado)

Se isso vale, U_ruim = conj(U_bom), logo |<x|U|0>|^2 e identico, logo o
oraculo de contagens tem poder EXATAMENTE zero -- nao aproximadamente.

Este script testa cada gate padrao contra essa condicao.
"""
import numpy as np

I2 = np.eye(2)
H  = np.array([[1,1],[1,-1]])/np.sqrt(2)
X  = np.array([[0,1],[1,0]], complex)
Y  = np.array([[0,-1j],[1j,0]])
Z  = np.diag([1,-1]).astype(complex)
S  = np.diag([1,1j]);           SDG = np.diag([1,-1j])
T  = np.diag([1,np.exp(1j*np.pi/4)]); TDG = np.diag([1,np.exp(-1j*np.pi/4)])
SX = 0.5*np.array([[1+1j,1-1j],[1-1j,1+1j]])

def rx(t): c,s=np.cos(t/2),np.sin(t/2); return np.array([[c,-1j*s],[-1j*s,c]])
def ry(t): c,s=np.cos(t/2),np.sin(t/2); return np.array([[c,-s],[s,c]], complex)
def rz(t): return np.diag([np.exp(-1j*t/2), np.exp(1j*t/2)])
def rz_bad(t): return np.diag([np.exp(1j*t/2), np.exp(-1j*t/2)])
def phase(t): return np.diag([1, np.exp(1j*t)])

def same_up_to_global(A, B, tol=1e-12):
    """A == e^{i phi} B ?  Fase global nao e observavel isoladamente."""
    idx = np.unravel_index(np.argmax(np.abs(B)), B.shape)
    if abs(B[idx]) < tol: return np.allclose(A, B, atol=tol)
    ph = A[idx]/B[idx]
    return abs(abs(ph)-1) < 1e-9 and np.allclose(A, ph*B, atol=tol)

print("A falha injetada e: inverter o sinal do Rz. Nada mais muda.")
print("Pergunta por gate: conj(g) == g (real), ou conj(g) == o que o backend defeituoso produz?\n")
print(f"{'gate':<10}{'conj = ele mesmo':<20}{'conj = versao defeituosa':<26}{'preserva a cegueira?'}")
print("-"*82)

fixos = [("I",I2),("H",H),("X",X),("Y",Y),("Z",Z),("S",S),("Sdg",SDG),
         ("T",T),("Tdg",TDG),("SX",SX)]
for nome, g in fixos:
    real  = np.allclose(np.conj(g), g)
    reala = same_up_to_global(np.conj(g), g)     # real a menos de fase global
    print(f"{nome:<10}{str(real):<20}{'—':<26}{'SIM' if real else ('so isolado' if reala else 'NAO')}")

ths = np.linspace(0.1, 3.0, 7)
for nome, f, fbad in [("Rx",rx,None),("Ry",ry,None),("Rz",rz,rz_bad),("P",phase,None)]:
    real = all(np.allclose(np.conj(f(t)), f(t)) for t in ths)
    if fbad:
        casa = all(np.allclose(np.conj(f(t)), fbad(t)) for t in ths)
        print(f"{nome:<10}{str(real):<20}{str(casa):<26}{'SIM' if casa else 'NAO'}")
    else:
        print(f"{nome:<10}{str(real):<20}{'—':<26}{'SIM' if real else 'NAO'}")
