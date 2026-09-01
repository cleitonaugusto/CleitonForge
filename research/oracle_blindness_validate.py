"""Valida a predicao analitica com circuitos aleatorios.

Predicao: sobre a familia {I,H,X,Z,Ry,Rz,CX,CZ,SWAP}, o oraculo de CONTAGENS
tem poder EXATAMENTE zero contra inversao de sinal do Rz. Sair da familia
(acrescentar S, T, SX, Rx ou P) deve quebrar a cegueira.
"""
import numpy as np
rng = np.random.default_rng(20260901)

I2=np.eye(2); H=np.array([[1,1],[1,-1]])/np.sqrt(2)
X=np.array([[0,1],[1,0]],complex); Y=np.array([[0,-1j],[1j,0]])
Z=np.diag([1,-1]).astype(complex); S=np.diag([1,1j]); T=np.diag([1,np.exp(1j*np.pi/4)])
SX=0.5*np.array([[1+1j,1-1j],[1-1j,1+1j]])
def rx(t): c,s=np.cos(t/2),np.sin(t/2); return np.array([[c,-1j*s],[-1j*s,c]])
def ry(t): c,s=np.cos(t/2),np.sin(t/2); return np.array([[c,-s],[s,c]],complex)
def rz(t,bad=False): e=1 if bad else -1; return np.diag([np.exp(1j*e*t/2),np.exp(-1j*e*t/2)])

N=4; DIM=2**N
def op1(g,q):
    m=np.array([[1]])
    for i in range(N): m=np.kron(m, g if i==q else I2)
    return m
def cx(c,t):
    M=np.zeros((DIM,DIM))
    for i in range(DIM):
        b=[(i>>(N-1-k))&1 for k in range(N)]
        if b[c]: b[t]^=1
        M[sum(v<<(N-1-k) for k,v in enumerate(b)), i]=1
    return M
def cz(a,b):
    M=np.eye(DIM,dtype=complex)
    for i in range(DIM):
        if ((i>>(N-1-a))&1) and ((i>>(N-1-b))&1): M[i,i]=-1
    return M

FAMILIAS = {
  "{H,X,Z,Ry,Rz,CX,CZ}  (predito CEGO)": ["H","X","Z","Ry","Rz","CX","CZ"],
  "  + Y                              ": ["H","X","Z","Ry","Rz","CX","CZ","Y"],
  "  + S                              ": ["H","X","Z","Ry","Rz","CX","CZ","S"],
  "  + T                              ": ["H","X","Z","Ry","Rz","CX","CZ","T"],
  "  + Rx                             ": ["H","X","Z","Ry","Rz","CX","CZ","Rx"],
  "  + SX                             ": ["H","X","Z","Ry","Rz","CX","CZ","SX"],
}

def build(gates, depth, bad):
    U=np.eye(DIM,dtype=complex)
    for _ in range(depth):
        g=gates[rng.integers(len(gates))]
        if g in ("CX","CZ"):
            a,b=rng.choice(N,2,replace=False)
            U=(cx(a,b) if g=="CX" else cz(a,b))@U
        else:
            q=int(rng.integers(N)); t=float(rng.uniform(-np.pi,np.pi))
            m={"H":H,"X":X,"Y":Y,"Z":Z,"S":S,"T":T,"SX":SX,
               "Rx":rx(t),"Ry":ry(t),"Rz":rz(t,bad)}[g]
            U=op1(m,q)@U
    return U

TRIALS, DEPTH = 400, 14
print(f"{TRIALS} circuitos aleatorios por familia, {N} qubits, profundidade {DEPTH}\n")
print(f"{'familia':<38}{'contagens diferem':<20}{'estado difere'}")
print("-"*76)
for nome, gates in FAMILIAS.items():
    dif_c = dif_s = 0
    for _ in range(TRIALS):
        st = rng.bit_generator.state
        a = build(gates, DEPTH, False)
        rng.bit_generator.state = st          # mesma sequencia, so o Rz muda
        b = build(gates, DEPTH, True)
        psi_a, psi_b = a[:,0], b[:,0]
        if not np.allclose(np.abs(psi_a)**2, np.abs(psi_b)**2, atol=1e-12): dif_c += 1
        if abs(abs(np.vdot(psi_a,psi_b))**2 - 1) > 1e-12: dif_s += 1
    print(f"{nome:<38}{dif_c/TRIALS*100:>6.1f}%{'':<13}{dif_s/TRIALS*100:>6.1f}%")
