import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.linalg import eigh

# __ Parameters _________________________________
L = 50         # length of the bar
E = 2.1e11      # Young's modulus
rho = 7850      # density
A = 0.01        # cross-sectional area

# __ Meshing ____________________________________
N = 500           # number of elements
nodes = np.linspace(0, L, N+1)
u = np.zeros(N+1)   # displacement matrix

# __ Stiffness Matrix ___________________________
K = np.zeros((N+1, N+1))
for i in range(N):
    K[i][i] += 1
    K[i+1][i] -= 1
    K[i][i+1] -= 1
    K[i+1][i+1] += 1

K *= E * A / (L/N)

# __ Mass Matrix ________________________________
M = np.zeros((N+1, N+1))
for i in range(N):
    M[i][i] += 2
    M[i+1][i] += 1
    M[i][i+1] += 1
    M[i+1][i+1] += 2
M *= rho * A * (L/N) / 6

# __ Boundary Conditions _________________________
K_red = K[1:, 1:]
M_red = M[1:, 1:]

# __ Eigenvalue Problem __________________________
modes = np.zeros((N+1, N+1))
freqs, modes[1:, 1:] = eigh(K_red, M_red)
freqs = np.sqrt(np.abs(freqs)) / (2 * np.pi)  # Convert to Hz
print(freqs[:10])  # Print first 10 natural frequencies

# __ Exporting Results ___________________________
df = pd.DataFrame({
    'Mode': np.arange(1, 11),
    'Frequency (Hz)': freqs[:10]
})
df.to_csv("axialvibrationsfrequency.csv", index = False)
# __ Plotting ____________________________________
plt.figure(figsize=(10, 6))
for i in range(4):  # Plot first 4 modes
    plt.plot(nodes, modes[:, i+1], label=f'Mode {i+1} (f={freqs[i]:.2f} Hz)')
plt.title('Axial Vibrations of a Bar')
plt.xlabel('Position along the bar (m)')
plt.ylabel('Displacement (arbitrary units)')
plt.legend()
plt.grid()
plt.show()
