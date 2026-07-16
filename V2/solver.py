import numpy as np
from scipy.linalg import eigh

# __ Internal Functions _______________
def ApplyBC(boundaryconditions:list):
    fixeddofs = []
    for bc in boundaryconditions:
        fixed = bc.fixeddof()
        fixeddofs += fixed
    fixeddofs = np.array(fixeddofs)

    return fixeddofs

def ApplyLoads(loading:list):
    loadlist = []
    for load in loading:
        if load.type == "pointload":
            if load.Fx != 0: loadlist.append((load.node.ID*3, load.Fx))
            if load.Fy != 0: loadlist.append((load.node.ID*3 + 1, load.Fy))
            if load.M != 0: loadlist.append((load.node.ID*3 + 2, load.M))
    
    return loadlist

def static(loading:list, boundaryconditions:list, K):
    N = K.shape[0] 
    alldofs = np.arange(0, N)
    fixeddofs = ApplyBC(boundaryconditions)
    freedofs = np.setdiff1d(alldofs, fixeddofs)

    K_red = K[np.ix_(freedofs, freedofs)]

    u = np.zeros(N)
    F = np.zeros_like(u)

    loadeddofs = ApplyLoads(loading)
    for load in loadeddofs: F[load[0]] = load[1]

    u[freedofs] = np.linalg.solve(K_red, F[np.ix_(freedofs)])

    return u, F

def modal(boundaryconditions:list, K, M):
    N = K.shape[0] 
    alldofs = np.arange(0, N)
    print(alldofs.shape)
    fixeddofs = ApplyBC(boundaryconditions)
    freedofs = np.setdiff1d(alldofs, fixeddofs)

    K_red = K[np.ix_(freedofs, freedofs)]
    M_red = M[np.ix_(freedofs, freedofs)]

    modes = np.zeros((N, N))
    freqs, modes[np.ix_(freedofs, freedofs)] = eigh(K_red, M_red)
    freqs = np.sqrt(np.abs(freqs)) / (2 * np.pi)

    return freqs, modes

# __ Still under progress _____________
def harmonic(loading:list, boundaryconditions:list, frequencies, K, M, C = False):
    N = K.shape[0] 
    alldofs = np.arange(0, N)
    print(alldofs.shape)
    fixeddofs = ApplyBC(boundaryconditions)
    freedofs = np.setdiff1d(alldofs, fixeddofs)

    F = u = np.zeros(N)
    loadeddofs = ApplyLoads(loading)
    for load in loadeddofs: F[load[0]] = load[1]
    
    if ~C:
        C = np.zeros_like(K)

    K_red = K[np.ix_(freedofs, freedofs)]
    M_red = M[np.ix_(freedofs, freedofs)]
    C_red = C[np.ix_(freedofs, freedofs)]
    F_red = F[np.ix_(freedofs)]

    
    U = []
    for f in frequencies: 
        w = 2 * np.pi * f
        D = K_red - (w**2)*M_red
        u = np.linalg.solve(D, F_red)
        U.append(u)
    
    return U
