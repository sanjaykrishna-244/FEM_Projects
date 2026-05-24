import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.linalg import eigh

# __ Parameters _________________________________
L = 50          # length of the bar
E = 2.1e11      # Young's modulus
rho = 7850      # density
A = 0.01        # cross-sectional area
I = 8.333e-6    # moment of inertia for a rectangular cross-section

# __ Meshing ____________________________________
N = 100                  # number of elements
nodes = np.linspace(0, L, N+1)
v = np.zeros(3*(N+1))   # displacement matrix
l = L/N

# __ Stiffness Matrix ___________________________
K = np.zeros((3*(N+1), 3*(N+1)))
a = 12
b = 6*l
c = 2*l**2
for i in range(0, 3*N, 3):
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
    k[np.ix_(axial, axial)] *= E*A/l
    k[np.ix_(bendg, bendg)] *= E*I/(l**3)
    K[i:i+6, i:i+6] += k


# __ Mass Matrix ________________________________
M = np.zeros((3*(N+1), 3*(N+1)))
for i in range(0, 3*N, 3):
    m = np.array([
                [2, 0, 0, 1, 0, 0],
                [0, 156, 22*l, 0, 54, -13*l],
                [0, 22*l, 4*l**2, 0, 13*l, -3*l**2],
                [1, 0, 0, 2, 0, 0],
                [0, 54, 13*l, 0, 156, -22*l],
                [0, -13*l, -3*l**2, 0, -22*l, 4*l**2]
                ])
    axial = [0, 3]
    bendg = np.setdiff1d(np.arange(0, 6), axial)
    m[np.ix_(axial, axial)] *= rho*A*l/6
    m[np.ix_(bendg, bendg)] *= rho*A*l/420
    M[i:i+6, i:i+6] += m

# __ Boundary Conditions _________________________
fixed = [0, 1, 2, 3*(N+1)-2, 3*(N+1)-1, 3*(N+1)]
free = np.setdiff1d(np.arange(3*(N+1)), fixed)
K_red = K[np.ix_(free, free)]
M_red = M[np.ix_(free, free)]

# __ Eigenvalue Problem __________________________
modes = np.zeros_like(K)
freqs, modes[np.ix_(free, free)] = eigh(K_red, M_red)
freqs = np.sqrt(np.abs(freqs)) / (2 * np.pi)  # Convert to Hz
#print(freqs[0: 25])  # Print first 25 natural frequencies

# __ Exporting Results ___________________________
df_freq = pd.DataFrame({
    'Mode': np.arange(1, 26),
    'Frequency (Hz)': freqs[0: 25]
})
df_space1 = pd.DataFrame([["Stiffness Matrix(K)"] + [""]*K.shape[1]])
df_stif = pd.DataFrame(K)
df_space2 = pd.DataFrame([[""]*K.shape[1], ["Masss Matrix(K)"] + [""]*K.shape[1]])
df_mass = pd.DataFrame(M)
df_mat = pd.concat([df_space1, df_stif, df_space2, df_mass], axis=0, ignore_index=True)
with pd.ExcelWriter("1Dbar/combinedaxialandbending.xlsx") as writer:
    df_freq.to_excel(writer, sheet_name='Frequencies', index=False)
    df_mat.to_excel(writer, sheet_name='Matrices', index=False, header=False)


    