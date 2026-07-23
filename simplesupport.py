import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import V2.geometry.crosssection as sec
import V2.geometry.definition as defin
import V2.material as mat

import V2.mesh as mesh
import V2.asssemble as ass

import V2.loads as lds
import V2.solver as sol

import V2.postprocess as pstprcs

# __ Material _________________________
steel = mat.Material()      # Using the default material

# __ Cross-section ____________________
section = sec.Rectangle(width=0.05, thickness=0.05)   # Using the default cross-section (0.01 X 0.01)

# __ Simple Supported Beam ____________
beam1 = defin.Beam_element(
    start=(0, 0),       # Start Point of Beam
    end = (0.5, 0),     # End Point of Beam
    material = steel,   # Material
    section = section   # Crossection
)

# __ Discretisation and Meshing _______
n = 40
frame = np.array([(beam1, n)])
nodes, elements, nodemap = mesh.OneDMeshing(frame)

# __ Assembling the matrices __________
K, M = ass.assemble(nodes, elements)

# __ Applying loads and boundary conditions __
simplesupportpoint1 = beam1.start
simplesupportpoint2 = beam1.end
loadpoint = (0.25, 0)

simplesupport1 = lds.PinSupport(nodemap[simplesupportpoint1])
simplesupport2 = lds.PinSupport(nodemap[simplesupportpoint2])
Load1 = lds.PointLoad(nodemap[loadpoint], Fy=-10000)

U, f = sol.static(loading=[Load1], boundaryconditions=[simplesupport1, simplesupport2], K=K)
#print(f"The deflection of simply supported beam at middle using FEM = {U[Load1.node.ID*3 + 1]}")


umax = (Load1.Fy * (beam1.length() ** 3)) / (48 * steel.E * section.Izz)
#print(f"The deflection of simply supported beam at middle with analytical solutions = {umax}")

err = (np.abs((U[Load1.node.ID*3 + 1] - umax) / umax)) * 100
#print(f"Percentage error = {err}")

results = pstprcs.EBBTpostprocess(U, elements, K, f)
strain = results["Strain"]
stress = results["Stress"]
reaction = results["Reactions"]
internalforce = results["Internal Forces"]

#print(stress[0].shape)

freqs, modes = sol.modal(boundaryconditions=[simplesupport1, simplesupport2], K = K, M = M)
#print(modes[np.arange(0, 123, 3), 7])
#print(len(freqs), modes.shape)
#print(nodes[0].coords)
XY = np.array([list(node.coords) for node in nodes])
#print(XY.shape, len(nodes))
x = XY[:, 0]
y = XY[:, 1]
i = int(input("Enter the mode number:"))
modei = modes[:, i]
u = modei[np.arange(0, 3*len(nodes), 3)]
v = modei[np.arange(1, 3*len(nodes), 3)]
X = x + (u )
Y = y + (v )

plt.figure()
plt.plot(x, y, color='k')
plt.plot(X, Y, color = (0, 1, 0))
plt.axis("equal")
plt.grid(True)
plt.show()

#print(modes.T @ M @ modes)

'''dfmodes = pd.DataFrame(modes)
with pd.ExcelWriter(path="sample.xlsx") as writer:
    dfmodes.to_excel(writer, sheet_name="Modeshapes")
'''
print(freqs[:12])
frequncies = np.arange(0, 20000, 10)
FRF = sol.harmonic_modalsuperposition(loading=[Load1], boundaryconditions=[simplesupport1, simplesupport2], frequencies=frequncies, K = K, M = M)

modei = FRF[:, 572]
u = modei[np.arange(0, 3*len(nodes), 3)]
v = modei[np.arange(1, 3*len(nodes), 3)]
X = x + (u )
Y = y + (v )

plt.figure()
plt.plot(x, y, color='k')
plt.plot(X, Y, color = (0, 1, 0))

plt.grid(True)
plt.show()

response22 = FRF[61, :]
plt.figure()
plt.plot(frequncies, response22)
plt.title("Modal Superposition")

plt.grid()
plt.show()

FRF = sol.harmonic_full(loading=[Load1], boundaryconditions=[simplesupport1, simplesupport2], frequencies=frequncies, K = K, M = M)

modei = FRF[:, 572]
u = modei[np.arange(0, 3*len(nodes), 3)]
v = modei[np.arange(1, 3*len(nodes), 3)]
X = x + (u )
Y = y + (v )

plt.figure()
plt.plot(x, y, color='k')
plt.plot(X, Y, color = (0, 1, 0))

plt.grid(True)
plt.show()

response22 = FRF[61, :]
plt.figure()
plt.plot(frequncies, response22)
plt.title("Full method")
plt.grid()
plt.show()
