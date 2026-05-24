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
v = np.zeros(2*(N+1))   # displacement matrix
l = L/N
# __ Stiffness Matrix ___________________________
K = np.zeros((2*(N+1), 2*(N+1)))
for i in range(0, 2*N, 2):
    k = np.array([[12, 6*l, -12, 6*l], 
                  [6*l, 4*l**2, -6*l, 2*l**2],
                  [-12, -6*l, 12, -6*l],
                  [6*l, 2*l**2, -6*l, 4*l**2]])
    K[i:i+4, i:i+4] += k
K *= E*I/(l**3)

# __ Mass Matrix ________________________________
M = np.zeros((2*(N+1), 2*(N+1)))
for i in range(0, 2*N, 2):
    m = np.array([[156, 22*l, 54, -13*l],
                  [22*l, 4*l**2, 13*l, -3*l**2],
                  [54, 13*l, 156, -22*l],
                  [-13*l, -3*l**2, -22*l, 4*l**2]])
    M[i:i+4, i:i+4] += m
M *= rho*A*l/420

# __ Boundary Conditions _________________________
fixed = [0, 2*N]
free = np.setdiff1d(np.arange(2*(N+1)), fixed)
print(free)  # Print free DOFs for verification
K_red = K[np.ix_(free, free)]
M_red = M[np.ix_(free, free)]


# __ Eigenvalue Problem __________________________
modes = np.zeros_like(K_red)
freqs, modes = eigh(K_red, M_red)
freqs = np.sqrt(np.abs(freqs)) / (2 * np.pi)  # Convert to Hz
print(freqs[:10])  # Print first 10 natural frequencies

# __ Normalisation of modes ______________________
'''eta = np.zeros((len(free), 4))
for i in range(2, 6):
    m = np.transpose(modes[:, i]) @ M_red @ modes[:, i]
    print(m)
    eta[:, i-2] = modes[:, i] / np.sqrt(m)
'''
# __ Exporting Results ___________________________
df = pd.DataFrame({
    'Mode': np.arange(1, 21),
    'Frequency (Hz)': freqs[:20]
})
df.to_csv("bendingvibrationsfrequency.csv", index = False)
ref_modes = np.zeros((2*(N+1), modes.shape[1]))
ref_modes[free, :] = modes
print(np.shape(ref_modes), np.size(nodes))

# __ Plotting _____________________________________
plt.figure(figsize=(10, 6))
for i in range(4):  # Plot first 4 modes
    plt.plot(nodes, ref_modes[0::2, i], label=f'Mode {i+1} (f={freqs[i]:.2f} Hz)')
plt.xlabel('Position along the bar')
plt.ylabel('Displacement')
plt.title('Bending Vibrations - First 4 Modes')
plt.legend()
plt.grid(True)
plt.savefig("bendingvibrationsmodes.png")
plt.show()
