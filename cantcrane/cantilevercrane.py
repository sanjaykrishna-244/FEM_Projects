import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.linalg import eigh
import utlities as ut

# __ Structure ________________________
"""
                | (3,4)
                |
                |
                | (3,0)
   (0,0)________|________ (6,0)
                \
                 \
                  \
                   \
            (5, -3) \ 

"""

# __ Properties _______________________
A = 0.01        # m2,    Area of cross-section
I = 8.3333e-6   # m4,    second moment of area
E = 210e+9      # Pa,    Young's Modulus
rho = 7850      # kg/m3, Density

# __ Geometric Properties _____________
L1 = np.sqrt((3-0)**2 + (0-0)**2); N1 = 20; dl1 = L1/N1
L2 = np.sqrt((6-3)**2 + (0-0)**2); N2 = 20; dl2 = L2/N2
L3 = np.sqrt((3-3)**2 + (4-0)**2); N3 = 25; dl3 = L3/N3
L4 = np.sqrt((5-3)**2 + (0-(-3))**2); N4 = 25; dl4 = L4/N4

# __ Discretization ___________________
# __ Member1 __:
x1 = np.linspace(0, L1, N1+1); x1_ = x1
y1 = np.zeros_like(x1); y1_ = y1

# __ Member2 __:
x2 = np.linspace(0, L2, N2+1) + L1; x2_ = x2[1:]
y2 = np.zeros_like(x2); y2_ = y2[1:]

# __ Member3 __:
y3 = np.linspace(0, L3, N3+1); y3_ = y3[1:]
x3 = np.zeros_like(y3) + 3; x3_ = x3[1:]

# __ Member4 __:
x4 = np.linspace(3, 5, N4+1); x4_ = x4[1:]
y4 = np.linspace(0, -3, N4+1); y4_ = y4[1:]

# __ all nodes __:
x = np.concat((x1_, x2_, x3_, x4_))
y = np.concat((y1_, y2_, y3_, y4_))
N = x.shape[0]

dofs = []
for i in range(N):
    j = 3*i
    dof = (j, j+1, j+2,)
    dofs.append(dof)
dofs = np.array(dofs)

# __ coordinate __:
df_coord = pd.DataFrame({
    "X":x,
    "Y":y,
    "DOF_u": dofs[:, 0],
    "DOF_v": dofs[:, 1],
    "DOF_t": dofs[:, 2]
})
"""
        DOFS                        RANGE
01 -> 21 = member1 nodes        => 0  -> 21 
21 -> 41 = member2 nodes        => 20 -> 41
21, 42 -> 66 = member3 nodes    => 20, 41 -> 66
21, 67 -> 91 = member4 nodes    => 20, 66 -> 91
"""

# __ Matrices _________________________
K = np.zeros((3*N, 3*N))
M = np.zeros((3*N, 3*N))


# __ Member1 __:
dof1 = np.arange(0, 21)
T1 = ut.transformation(0)
for i in range(len(dof1) - 1):
    k = ut.stiffness(T1, dl1, E, A, I)
    m = ut.mass(T1, dl1, A, rho)
    ut.matrixform(k, K, i, dof1)
    ut.matrixform(m, M, i, dof1)

# __ Member2 __:
dof2 = np.arange(20, 41)
T2 = ut.transformation(0)
for i in range(len(dof2) - 1):
    k = ut.stiffness(T2, dl2, E, A, I)
    m = ut.mass(T2, dl2, A, rho)
    ut.matrixform(k, K, i, dof2)
    ut.matrixform(m, M, i, dof2)

# __ Member3 __:
dof3 = np.concat((np.array([20,]), np.arange(41, 66)))
T3 = ut.transformation(90)
for i in range(len(dof3) - 1):
    k = ut.stiffness(T3, dl3, E, A, I)
    m = ut.mass(T3, dl3, A, rho)
    ut.matrixform(k, K, i, dof3)
    ut.matrixform(m, M, i, dof3)

# __ Member4 __:
ang4 = np.rad2deg(np.atan2((-3 - 0), (5 - 3)))
print(ang4)
dof4 = np.concat((np.array([20,]), np.arange(66, 91)))
T4 = ut.transformation(ang4)
for i in range(len(dof4) - 1):
    k = ut.stiffness(T4, dl4, E, A, I)
    m = ut.mass(T4, dl4, A, rho)
    ut.matrixform(k, K, i, dof4)
    ut.matrixform(m, M, i, dof4)

dfname = pd.DataFrame([["Stiffness_Matrix"] + [""]*(K.shape[1]-1)])
df_k = pd.DataFrame(K)
dfk = pd.concat([dfname, df_k], ignore_index=True)

dfname = pd.DataFrame([["Mass Matrix"]+[""]*(M.shape[1]-1)])
df_m = pd.DataFrame(M)
dfm = pd.concat([dfname, df_m], ignore_index=True)

# __ Boundary Conditions ______________
fixed = np.array([])
free = np.setdiff1d(np.arange(0, 3*N), fixed)
K_red = K[np.ix_(free, free)]
M_red = M[np.ix_(free, free)]

# __ Eigen Value Problem ______________
modes = np.zeros_like(K)
freqs, modes[np.ix_(free, free)] = eigh(K_red, M_red)
freqs = np.sqrt(np.abs(freqs)) / 2*np.pi
dff = pd.DataFrame({
    'Mode' : np.arange(1,31),
    'Frequency' : freqs[:30]
})

# __ Plotting _________________________
i = 0
plt.figure()
plt.plot(x1, y1, 'b'); plt.plot(x2, y2, 'b'); plt.plot(x3, y3, 'b'); plt.plot(x4, y4, 'b'); 
plt.plot(x1 + 5* modes[:, i][np.ix_(dof1*3)], y1 + 5* modes[:, i][np.ix_(dof1 * 3 + 1)], 'g')
plt.plot(x2 + 5* modes[:, i][np.ix_(dof2*3)], y2 + 5* modes[:, i][np.ix_(dof2 * 3 + 1)], 'g')
plt.plot(x3 + 5* modes[:, i][np.ix_(dof3*3)], y3 + 5* modes[:, i][np.ix_(dof3 * 3 + 1)], 'g')
plt.plot(x4 + 5* modes[:, i][np.ix_(dof4*3)], y4 + 5* modes[:, i][np.ix_(dof4 * 3 + 1)], 'g')
plt.title(f"Cantilver crane {i-2}th mode shape")
plt.grid()
plt.show()

# __ Exporting results ________________
with pd.ExcelWriter("scantcrane/cranecantilver.xlsx") as writer:
    df_coord.to_excel(writer, sheet_name="Coordinates of Frame")
    dfk.to_excel(writer, sheet_name="Stiffness", header = False, index = False)
    dfm.to_excel(writer, sheet_name="Mass", header = False, index = False)
    dff.to_excel(writer, sheet_name="Frequencies")
