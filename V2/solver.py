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
    # Loading List contains all the loads along with the loactions as list of loading objects
    # boundary conditions is the list of all the boundary condition objects(fixed, pin, roller)
    # K is the assembled stiffness matrix without any dof reductions
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
    # M is the assembled mass matrix without any dof reductions
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
def harmonic_full(loading:list, boundaryconditions:list, frequencies, K, M, C = None):
    # frequencies are the sweeping set of frequencies for FRF calculation
    N = K.shape[0] 
    alldofs = np.arange(0, N)
    print(alldofs.shape)
    fixeddofs = ApplyBC(boundaryconditions)
    freedofs = np.setdiff1d(alldofs, fixeddofs)

    F = u = np.zeros(N)
    loadeddofs = ApplyLoads(loading)
    for load in loadeddofs: F[load[0]] = load[1]
    
    if C is None:
        C = np.zeros_like(K)

    K_red = K[np.ix_(freedofs, freedofs)]
    M_red = M[np.ix_(freedofs, freedofs)]
    C_red = C[np.ix_(freedofs, freedofs)]
    F_red = F[np.ix_(freedofs)]

    
    responses = []
    for f in frequencies: 
        w = 2 * np.pi * f
        D = K_red - (w**2)*M_red
        u = np.linalg.solve(D, F_red)
        responses.append(u)

    responses = np.array(responses).T
    U = np.zeros((N, len(frequencies)))
    U[freedofs, :] = responses
    
    return U

def harmonic_modalsuperposition(loading:list, boundaryconditions:list, frequencies, K, M, C = None, modeshapes = None, retainedmodes = 10):
    # Modesahpes matrix of column vectors of modes shapes that are direct solution of modal analysis
    # naturalfreq: set of retained natural frequencies whose modeshape are used to calculate the FRF
    N = K.shape[0] 
    alldofs = np.arange(0, N)
    print(alldofs.shape)
    fixeddofs = ApplyBC(boundaryconditions)
    freedofs = np.setdiff1d(alldofs, fixeddofs)

    F = u = np.zeros(N)
    loadeddofs = ApplyLoads(loading)
    for load in loadeddofs: F[load[0]] = load[1]
    
    if C is None:
        C = np.zeros_like(K)

    if modeshapes is None:
        _, modeshapes = modal(boundaryconditions, K, M)
        modeshapes = modeshapes[np.ix_(freedofs, freedofs)]

    retainedmodes = min(retainedmodes, len(freedofs))

    K_red = K[np.ix_(freedofs, freedofs)]
    M_red = M[np.ix_(freedofs, freedofs)]
    C_red = C[np.ix_(freedofs, freedofs)]
    F_red = F[np.ix_(freedofs)]
    phi = modeshapes[:, :retainedmodes]

    Km = phi.T @ K_red @ phi; Km = np.round(Km, 8)
    Mm = phi.T @ M_red @ phi; Mm = np.round(Mm, 8)
    Cm = phi.T @ C_red @ phi; Cm = np.round(Cm, 8)
    Fm = phi.T @ F_red

    k = np.diag(Km)
    m = np.diag(Mm)
    c = np.diag(Cm)
    responses = []

    for i, f in enumerate(frequencies):
        w = 2 * np.pi * f
        D = k - (w ** 2) * m
        q = Fm / D
        u = phi @ q
        responses.append(u)

    responses = np.array(responses).T
    U = np.zeros((N, len(frequencies)))
    U[freedofs, :] = responses

    return U
    





