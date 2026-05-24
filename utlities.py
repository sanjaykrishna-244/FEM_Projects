import numpy as np

# __ Matrix Formations ________________
# __ Tansformation __:
def transformation(theta):
    theta = np.deg2rad(theta)
    c = np.cos(theta) if abs(np.cos(theta)) > 1e-10 else 0
    s = np.sin(theta) if abs(np.sin(theta)) > 1e-10 else 0
    T = np.array([[c, s, 0, 0, 0, 0],
                  [-s, c, 0, 0, 0, 0],
                  [0, 0, 1, 0, 0, 0],
                  [0, 0, 0, c, s, 0],
                  [0, 0, 0, -s, c, 0],
                  [0, 0, 0, 0, 0, 1]])
    return T
# __ Stiffness __:
def stiffness(T, dl, E, A, I):
    a = 12; b = 6 * dl; c = 2 * dl ** 2
    k = np.array([
                [1, 0, 0, -1, 0, 0],
                [0, a, b, 0, -a, b], 
                [0, b, 2*c, 0, -b, c], 
                [-1, 0, 0, 1, 0, 0], 
                [0, -a, -b, 0, a, -b], 
                [0, b, c, 0, -b, 2*c]
                ])
    axial = [0, 3]
    bendg = np.setdiff1d(np.arange(0, 6), axial)
    k[np.ix_(axial, axial)] *= E*A/dl
    k[np.ix_(bendg, bendg)] *= E*I/(dl**3)
    k = np.transpose(T) @ k @ T

    return k
# __ Mass __:
def mass(T, dl, A, rho):
    m = np.array([
        [140, 0, 0, 70, 0, 0],
        [0, 156, 22*dl, 0, 54, -13*dl],
        [0, 22*dl, 4*dl**2, 0, 13*dl, -3*dl**2],
        [70, 0, 0, 140, 0, 0],
        [0, 54, 13*dl, 0, 156, -22*dl],
        [0, -13*dl, -3*dl**2, 0, -22*dl, 4*dl**2]]
        )
    m *= rho * A * dl / 420
    m = T.T @ m @ T

    return m
# __ Final matrix addition __:
def matrixform(m, M, i, D):
    i1 = 3*D[i]
    i2 = 3*D[i+1]
    dof = np.array([i1, i1+1, i1+2, i2, i2+1, i2+2])
    M[np.ix_(dof, dof)] += m
