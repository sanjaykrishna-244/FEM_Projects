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

# __ Material _________________________
steel = mat.Material()

# __ Cross-section ____________________
section = sec.Rectangle(width=0.05, thickness=0.05)

# __ Combination of Beams _____________
beam1 = defin.StraightBeam_element(
    start = (0, 0),
    end = (1, 0),
    material = steel,
    section = section
); N1 = 100
beam2 = defin.StraightBeam_element(
    start = (1, 0),
    end = (2, 0),
    material = steel,
    section = section
); N2 = 100
beam3 = defin.StraightBeam_element(
    start = (0, 1),
    end = (1, 1),
    material = steel,
    section = section
); N3 = 100
beam4 = defin.StraightBeam_element(
    start = (1, 1),
    end = (2, 1),
    material = steel,
    section = section
); N4 = 100
beam5 = defin.StraightBeam_element(
    start = (0, 0),
    end = (0, 1),
    material = steel,
    section = section
); N5 = 100
beam6 = defin.StraightBeam_element(
    start = (1.5, 0),
    end = (1.5, 1),
    material = steel,
    section = section
); N6 = 100
beam7 = defin.StraightBeam_element(
    start = (2, 0),
    end = (2, 1),
    material = steel,
    section = section
); N7 = 100

# __ Discretisation and Meshing _______
frame = np.array([(beam1, N1), (beam2, N2), (beam3, N3), (beam4, N4), (beam5, N5), (beam6, N6), (beam7, N7)])
nodes, elements, nodesmap = mesh.OneDMeshing(frame)

# __ Assembling the matrices __________
K, M = ass.assemble(nodes, elements)

plt.figure()
for element in elements:
    x = [element.node1.coords[0], element.node2.coords[0]]
    y = [element.node1.coords[1], element.node2.coords[1]]
    plt.plot(x, y, color = (1, 0, 1))
plt.axis('equal')
plt.grid()
plt.title("Frame Without Loading")
plt.xlabel("in (m)")
plt.ylabel("in (m)")
plt.savefig('deckframe/nodeformation.png')
plt.show()

# __ Applying loads and boundary conditions __
loadpoint1 = (1.5, 0)
simplesupportpoint1 = beam1.start
simplesupportpoint2 = beam2.end

simplesupport1 = lds.PinSupport(nodesmap[simplesupportpoint1])
simplesupport2 = lds.PinSupport(nodesmap[simplesupportpoint2])
load1 = lds.PointLoad(nodesmap[loadpoint1], Fy=-100000)

loadlist = [load1]
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
plt.title("Static Deformation under static load of -1000000 vertical load")
plt.xlabel("in (m)")
plt.ylabel("in (m)")
plt.savefig('deckframe/staticdeformation.png')
plt.show()


modalresults = sol.modal(boundaryconditions = boundaryconditionslist, K = K, M = M)
modeshapes = modalresults["Modeshapes"]
nat_freqs = modalresults["Natural Frequencies"]
print((nat_freqs[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]]))
plt.figure()
i = 2
modeshape1 = modeshapes[:, i]
for element in elements:
    dofs = element.dof()
    u = modeshape1[dofs[[0, 3]]]
    v = modeshape1[dofs[[1, 4]]]
    x = [element.node1.coords[0], element.node2.coords[0]] + u
    y = [element.node1.coords[1], element.node2.coords[1]] + v
    plt.plot(x, y, color = (1, 0, 0.5))
plt.grid()
plt.title(f"Modeshape of {i+1}th natural frequncy \n f = {nat_freqs[i]}")
plt.xlabel("in (m)")
plt.ylabel("in (m)")
plt.savefig('deckframe/modeshape3.png')
plt.show()
'''
freqs = np.concat(
    (np.arange(0, (nat_freqs[0] // 10) *10, 10),
    np.arange((nat_freqs[0] // 10) * 10, (nat_freqs[0] // 10) * 10 + 10, 0.1),
    np.arange((nat_freqs[0] // 10) * 10 + 10, (nat_freqs[1] // 10) * 10, 10),
    np.arange((nat_freqs[1] // 10) * 10, (nat_freqs[1] // 10) * 10 + 10, 0.1),
    np.arange((nat_freqs[1] // 10) * 10 + 10, (nat_freqs[2] // 10) * 10, 10),
    np.arange((nat_freqs[2] // 10) * 10, (nat_freqs[2] // 10) * 10 + 10, 0.1),
    np.arange((nat_freqs[2] // 10) * 10 + 10, 100, 10),)
)
harmonicresults = sol.harmonic_full(loadlist, boundaryconditionslist, frequencies= freqs, K = K, M = M)
responses = harmonicresults["Deformations"]
midspanVresponse = responses[nodesmap[loadpoint1].ID*3 + 1, :]
plt.figure()
plt.semilogy(freqs, np.abs(midspanVresponse), color = (0, 0, 1))
plt.title("FRF of simply supported beam\n at midspan under -1000000N harmonic load")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Deformation (m)")
plt.grid('both')
plt.savefig('deckframe/FRF.png')
plt.show()
'''