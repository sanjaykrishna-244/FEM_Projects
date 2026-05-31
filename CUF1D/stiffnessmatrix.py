import numpy as np
import sympy as sp

x, y, z = sp.symbols('x y z')
L = 1.0
E = 210e+9
n = 0.33
G = E / 2 * (1 + n)
l = 2 * G / (2 - n)
F1 = 1; F2 = y; F3 = z
F = np.array([[F1, F2, F3,  0,  0,  0,  0,  0,  0],
              [ 0,  0,  0, F1, F2, F3,  0,  0,  0],
              [ 0,  0,  0,  0,  0,  0, F1, F2, F3]])
N1 = 0.5 - 0.5*x; N2 = 0.5 + x

def Bi(F, Ni):
    N = np.eye(3) * Ni
    NF = N @ F
    B = np.array([
        np.array(sp.diff(sp.Matrix(NF[0]), x)).reshape(1, -1)[0], 
        np.array(sp.diff(sp.Matrix(NF[1]), y)).reshape(1, -1)[0], 
        np.array(sp.diff(sp.Matrix(NF[2]), z)).reshape(1, -1)[0], 
        np.array(sp.diff(sp.Matrix(NF[0]), y)).reshape(1, -1)[0] + np.array(sp.diff(sp.Matrix(NF[1]), x)).reshape(1, -1)[0], 
        np.array(sp.diff(sp.Matrix(NF[1]), z)).reshape(1, -1)[0] + np.array(sp.diff(sp.Matrix(NF[2]), y)).reshape(1, -1)[0], 
        np.array(sp.diff(sp.Matrix(NF[2]), x)).reshape(1, -1)[0] + np.array(sp.diff(sp.Matrix(NF[0]), z)).reshape(1, -1)[0], 
    ]
    )
    return B

C = np.array(
    [[l + G, G, G, 0, 0, 0],
     [G, l + G, G, 0, 0, 0],
     [G, G, l + G, 0, 0, 0],
     [0,   0,   0, G, 0, 0],
     [0,   0,   0, 0, G, 0],
     [0,   0,   0, 0, 0, G]]
)
B1 = Bi(F, N1)
B2 = Bi(F, N2)
B = np.array([B1, B2])
def Ke(B, C):
    ke = np.zeros((18, 18))
    print(type(ke[0][0]))
    for i in range(len(B)):
        for j in range(len(B)):
            BI = B[i]; BJ = B[j]
            BtCB = BJ.T @ C @ BI
            kij = sp.integrate(sp.integrate(sp.integrate(sp.Matrix(BtCB), (x, -0.5, 0.5)), (y, -0.05, 0.05)), (z, -0.05, 0.05))
            kij = np.array(kij)
            print(type(kij[0][0]))
            dofsi = np.arange(i*9, (i+1)*9)
            dofsj = np.arange(j*9, (j+1)*9)
            print(dofsi, dofsj)
            ke[np.ix_(dofsi, dofsj)] += kij
    print(ke.shape)
Ke(B, C)