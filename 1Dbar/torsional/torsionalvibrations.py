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
P = 0.33        # Poisson's ratio
G = E / (2 * (1 + P))   # Shear modulus, assuming Poisson's ratio of 0.33
J_p = (0.1**4) / 6      # Polar moment of inertia for a rectangular cross-section
J_t = 0.1406 * (0.1**4)  # Torsional constant for a rectangular cross-section

# __ Meshing ____________________________________
N = 500                  # number of elements
nodes = np.linspace(0, L, N+1)
t = np.zeros(N+1)   # displacement matrix
l = L/N

# __ Stiffness Matrix ___________________________
K = np.zeros(((N+1), (N+1)))
for i in range(N):
    k = np.array([[1, -1],
                 [-1, 1]])
    K[i:i+2, i:i+2] += k

K *= G*J_t/l

# __ Mass Matrix ________________________________
M = np.zeros(((N+1), (N+1)))
for i in range(N):
    m = np.array([[2, 1],
                 [1, 2]])
    M[i:i+2, i:i+2] += m
M *= rho * J_p * l / 6

# __ Boundary Conditions _________________________
fixed = [0,]
free = np.setdiff1d(np.arange(N+1), fixed)
#print(free)
K_red = K[np.ix_(free, free)]
M_red = M[np.ix_(free, free)]

# __ Eigenvalue Problem __________________________
modes = np.zeros((N+1, N+1))
freqs, modes[free[:, None], free] = eigh(K_red, M_red)
freqs = np.sqrt(np.abs(freqs)) / (2 * np.pi)  # Convert to Hz
print(freqs[:10])  # Print first 10 natural frequencies

# __ Exporting Results ___________________________
df = pd.DataFrame({
    'Mode': np.arange(1, 11),
    'Frequency (Hz)': freqs[:10]
})
df.to_csv("torsionalvibrationsfrequency.csv", index = False)

# __ Plotting ____________________________________
plt.figure(figsize=(10, 6))
for i in range(4):  # Plot first 4 modes
    plt.plot(nodes, modes[:, i+1], label=f'Mode {i+1} (f={freqs[i]:.2f} Hz)')
plt.title("Torsional Frequencies of a Bar with square cross-section")
plt.xlabel("Position along the bar(m)")
plt.ylabel("Torsional Displacement")
plt.legend()
plt.grid()
plt.show()
