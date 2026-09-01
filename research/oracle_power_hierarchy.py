"""Hierarquia de poder: contagens < estado < operador.

O oraculo de ESTADO ve uma coluna do operador (a partir de |0..0>).
O de OPERADOR ve todas. A diferenca deveria aparecer nos gates que agem como
identidade sobre o estado alcancado -- exatamente os 35% que a campanha do
Cirq perdia.
"""
import numpy as np
exec(open("validate.py").read().split("FAMILIAS =")[0])

FAM = {
  "{H,X,Z,Ry,Rz,CX,CZ}": ["H","X","Z","Ry","Rz","CX","CZ"],
  "+ T":                 ["H","X","Z","Ry","Rz","CX","CZ","T"],
  "+ Rx":                ["H","X","Z","Ry","Rz","CX","CZ","Rx"],
}
def build(gates, depth, bad):
    U=np.eye(DIM,dtype=complex)
    for _ in range(depth):
        g=gates[rng.integers(len(gates))]
        if g in ("CX","CZ"):
            a,b=rng.choice(N,2,replace=False); U=(cx(a,b) if g=="CX" else cz(a,b))@U
        else:
            q=int(rng.integers(N)); t=float(rng.uniform(-np.pi,np.pi))
            m={"H":H,"X":X,"Y":Y,"Z":Z,"S":S,"T":T,"SX":SX,
               "Rx":rx(t),"Ry":ry(t),"Rz":rz(t,bad)}[g]
            U=op1(m,q)@U
    return U

def op_equiv(A,B,tol=1e-10):
    """Operadores iguais a menos de fase global?"""
    i=np.unravel_index(np.argmax(np.abs(B)),B.shape)
    if abs(B[i])<tol: return np.allclose(A,B,atol=tol)
    return np.allclose(A, (A[i]/B[i])*B, atol=tol)

T_,D=600,14
print(f"{T_} circuitos por familia, {N} qubits, profundidade {D}")
print("poder = fracao de circuitos em que o oraculo DETECTA a inversao de sinal\n")
print(f"{'familia':<24}{'contagens':>11}{'estado':>10}{'operador':>11}")
print("-"*56)
res={}
for nome,gates in FAM.items():
    c=s=o=0
    for _ in range(T_):
        st=rng.bit_generator.state
        A=build(gates,D,False); rng.bit_generator.state=st; B=build(gates,D,True)
        pa,pb=A[:,0],B[:,0]
        if not np.allclose(np.abs(pa)**2,np.abs(pb)**2,atol=1e-12): c+=1
        if abs(abs(np.vdot(pa,pb))**2-1)>1e-12: s+=1
        if not op_equiv(A,B): o+=1
    res[nome]=(c/T_,s/T_,o/T_)
    print(f"{nome:<24}{c/T_*100:>10.1f}%{s/T_*100:>9.1f}%{o/T_*100:>10.1f}%")

print("\nIC 95% (Wilson) para a familia cega, oraculo de contagens:")
import math
n=T_; p=res["{H,X,Z,Ry,Rz,CX,CZ}"][0]; z=1.96
den=1+z*z/n; c1=(p+z*z/(2*n))/den
half=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
print(f"  {p*100:.1f}%  [{max(0,(c1-half))*100:.2f}%, {(c1+half)*100:.2f}%]")
