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
steel = mat.Material(E = 210e9, nu = 0.3, rho = 7850)
alum = mat.Material(E = 70e9, nu = 0.33, rho = 2700)

# __ Cross-section ____________________
sec1 = sec.Rectangle(width = 2.5e-3, thickness = 2e-3)
sec2 = sec.Rectangle(width = 2.5e-3, thickness = 1e-3)

# __ Combination of Beams _____________
beam1 = defin.StraightBeam_element(
    start = (0, 0),
    end = (0.6, 0),
    material = steel,
    section = sec1
); N1 = 100
beam2 = defin.StraightBeam_element(
    start = (0, 0.1),
    end = (0.6, 0.1),
    material = steel,
    section = sec1
); N2 = 100


# __ Meshing ___________________________
frame = [(beam1, N1), (beam2, N2)]
nodes, elements, nodesmap = mesh.OneDMeshing(frame)

plt.figure()
for element in elements:
    x = [element.node1.coords[0], element.node2.coords[0]]
    y = [element.node1.coords[1], element.node2.coords[1]]
    plt.plot(x, y, color = (1, 0, 1))
plt.axis('equal')
plt.grid()
plt.title("System Without Loading")
plt.xlabel("in (m)")
plt.ylabel("in (m)")
plt.savefig('btpmodel/no_deformation.png')
plt.show()