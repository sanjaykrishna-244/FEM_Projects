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
n = 100
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
Load1 = lds.PointLoad(nodemap[loadpoint], Fy=-1000000)

loadlist = [Load1]
boundaryconditionslist = [simplesupport1, simplesupport2]

staticresults = sol.static(loading=loadlist, boundaryconditions=boundaryconditionslist, K = K)
staticdeformation = staticresults["Deformation"]
plt.figure()
for element in elements:
    dofs = element.dof()
    u = staticdeformation[dofs[[0, 3]]]
    v = staticdeformation[dofs[[1, 4]]]
    x = [element.node1.coords[0], element.node2.coords[0]] + u
    y = [element.node1.coords[1], element.node2.coords[1]] + v
    plt.plot(x, y, color = (1, 0, 1))
plt.axis('equal')
plt.grid()
plt.title("Static Deformation under static load of -1000000 vertical load at midspan")
plt.xlabel("in (m)")
plt.ylabel("in (m)")
plt.savefig('staticdeformation.png')
plt.show()

modalresults = sol.modal(boundaryconditions = boundaryconditionslist, K = K, M = M)
modeshapes = modalresults["Modeshapes"]
nat_freqs = modalresults["Natural Frequencies"]
print((nat_freqs[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]] // 10 )* 10)
plt.figure()
modeshape1 = modeshapes[:, 2]
for element in elements:
    dofs = element.dof()
    u = modeshape1[dofs[[0, 3]]]
    v = modeshape1[dofs[[1, 4]]]
    x = [element.node1.coords[0], element.node2.coords[0]] + u
    y = [element.node1.coords[1], element.node2.coords[1]] + v
    plt.plot(x, y, color = (1, 0, 0.5))
plt.grid()
plt.title(f"Modeshape of 3rd natural frequncy \n f = {nat_freqs[2]}")
plt.xlabel("in (m)")
plt.ylabel("in (m)")
plt.savefig('modeshape3.png')
plt.show()

freqs = np.concat(
    (np.arange(0, (nat_freqs[0] // 10) *10, 10),
    np.arange((nat_freqs[0] // 10) * 10, (nat_freqs[0] // 10) * 10 + 10, 0.1),
    np.arange((nat_freqs[0] // 10) * 10 + 10, (nat_freqs[1] // 10) * 10, 10),
    np.arange((nat_freqs[1] // 10) * 10, (nat_freqs[1] // 10) * 10 + 10, 0.1),
    np.arange((nat_freqs[1] // 10) * 10 + 10, (nat_freqs[2] // 10) * 10, 10),
    np.arange((nat_freqs[2] // 10) * 10, (nat_freqs[2] // 10) * 10 + 10, 0.1),
    np.arange((nat_freqs[2] // 10) * 10 + 10, 5000, 10),)
)
harmonicresults = sol.harmonic_full(loadlist, boundaryconditionslist, frequencies= freqs, K = K, M = M)
responses = harmonicresults["Deformations"]
midspanVresponse = responses[nodemap[loadpoint].ID*3 + 1, :]
plt.figure()
plt.semilogy(freqs, np.abs(midspanVresponse), color = (0, 0, 1))
plt.title("FRF of simply supported beam\n at midspan under -1000000N harmonic load")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Deformation (m)")
plt.grid('both')
plt.savefig('FRF.png')
plt.show()