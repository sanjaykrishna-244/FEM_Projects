import numpy as np
import pandas as pd
# __ 1D 2 NODES ELEMENT STIFFNESS MATRIX LE Class __
#       Each element has 2 nodes (i = 1, 2)
#       Each Node is a cross-section with 4 nodes
#           u(x, y, z) = Ni(x) . Ft(y, z) . uit; Ni and Ft are shape functions
#           u and each uit has three components one each coordinate direction
#       Ni interpolates in x direction with 2 node cross-section
#       Ft interpolates the displacement field in the cross-section (y, z) with 4 nodes
#       
#       Now, (y, z) cross-section of each node has been normalised to (ξ, η)
#       Thus, the shape functions are defined in terms of ξ and η, and the element length is L
#               Lt(ξ, η) = (1 + ξt*ξ) * (1 + ηt*η) / 4; t = 1, 2, 3, 4; ηt = -1, 1 & v ξt = -1, 1
#               Ni(x) = (1/2 + (-1)**i * x / L); i = 1, 2                
#                   {-1 < η, ξ < 1; -L/2 < x < L/2}


# __ Material Properties _________________
E = 210e+9
n = 0.33
G = E / (2 * (1 + n))
l = 2 * G / (1 - 2 * n)

# __ Element Geometry _________________
L = 1

# __ Stiffness Matrix _________________
C = np.array([
    [l + 2 * G, l, l, 0, 0, 0],
    [l, l + 2 * G, l, 0, 0, 0],
    [l, l, l + 2 * G, 0, 0, 0],
    [    0,   0,   0, G, 0, 0],
    [    0,   0,   0, 0, G, 0],
    [    0,   0,   0, 0, 0, G],
])

# __ Useful Functions _________________
def FN(i, coord):
    if i == 1:
        N = 0.5 - (coord[0] / L)
    else: 
        N = 0.5 + (coord[0] / L)

    NF = 0.25 * np.array([
        N * (1 - coord[1]) * (1 + coord[2]),
        N * (1 + coord[1]) * (1 + coord[2]),
        N * (1 + coord[1]) * (1 - coord[2]),
        N * (1 - coord[1]) * (1 - coord[2]),
    ])

    return NF

def dFN_dx(i, coord):
    if i == 1:
        N = - 1 / L
    else: 
        N = + 1 / L
    
    NF = 0.25 *  np.array([
        N * (1 - coord[1]) * (1 + coord[2]),
        N * (1 + coord[1]) * (1 + coord[2]),
        N * (1 + coord[1]) * (1 - coord[2]),
        N * (1 - coord[1]) * (1 - coord[2]),
    ])

    return NF

def dFN_dzeta(i, coord):
    if i == 1:
        N = 0.5 - (coord[0] / L)
    else: 
        N = 0.5 + (coord[0] / L)
    
    NF = 0.25 * np.array([
        N * -1 * (1 + coord[2]),
        N * 1 * (1 + coord[2]),
        N * 1 * (1 - coord[2]),
        N * -1 * (1 - coord[2]),
    ])

    return NF

def dFN_deta(i, coord):
    if i == 1:
        N = 0.5 - (coord[0] / L)
    else: 
        N = 0.5 + (coord[0] / L)
    
    NF = 0.25 * np.array([
        N * (1 - coord[1]) * 1,
        N * (1 + coord[1]) * 1,
        N * (1 + coord[1]) * -1,
        N * (1 - coord[1]) * -1,
    ])

    return NF

def b(i, coord):
    zero = np.array([0, 0, 0, 0])
    B = np.array([
        np.concatenate([dFN_dx(i, coord), zero, zero]),
        np.concatenate([zero, dFN_dzeta(i, coord), zero]),
        np.concatenate([zero, zero, dFN_deta(i, coord)]),
        np.concatenate([zero, dFN_deta(i, coord), dFN_dzeta(i, coord)]),
        np.concatenate([dFN_deta(i, coord), zero, dFN_dx(i, coord)]),
        np.concatenate([dFN_dzeta(i, coord), dFN_dx(i, coord), zero]),
    ])

    return B

def ke(coord, C):
    Ke = np.zeros((24, 24))
    for i in range(2):
        for j in range(2):
            Bi = b(i + 1, coord)
            Bj = b(j + 1, coord)
            n = Bi.shape[1]
            kij = Bi.T @ C @ Bj
            dofi = np.arange(i * n, (i + 1) * n)
            dofj = np.arange(j * n, (j + 1) * n)
            Ke[np.ix_(dofi, dofj)] += kij

    return Ke


# __ Elemental Stiffness Matrix _______
δ = 1 / np.sqrt(3)
gp = np.array([
    [-δ*L/2, -δ,  δ],
    [-δ*L/2,  δ,  δ],
    [-δ*L/2,  δ, -δ],
    [-δ*L/2, -δ, -δ],
    [ δ*L/2, -δ,  δ],
    [ δ*L/2,  δ,  δ],
    [ δ*L/2,  δ, -δ],
    [ δ*L/2, -δ, -δ],
])
wt = np.array([1, 1, 1, 1, 1, 1, 1, 1])

intgrl = np.array([ke(coord, C) * wt[i] for i, coord in enumerate(gp)])
intgrl = np.sum(intgrl, axis=0)
df = pd.DataFrame(intgrl)
df.to_excel("CUF1D/LEclassstiffness.xlsx", index=False, header=False)
print(np.allclose(intgrl, intgrl.T))
print(np.linalg.eigvalsh(intgrl))
print(np.max(np.abs(intgrl - intgrl.T)))

'''for i, coord in enumerate(gp):
    B = b(1, coord)
    eig = np.linalg.eigvalsh(B.T @ C @ B)
    print(f"{i + 1}: {eig}", end="\n")
print("##########################################")
for i, coord in enumerate(gp):
    K = ke(coord, C)
    eig = np.linalg.eigvalsh(K)
    print(f"{i + 1}: {eig}", end="\n")'''

