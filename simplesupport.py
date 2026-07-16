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
section = sec.Rectangle()   # Using the default cross-section (0.01 X 0.01)

# __ Simple Supported Beam ____________
beam1 = defin.Beam_element(
    start=(0, 0),       # Start Point of Beam
    end = (0.5, 0),     # End Point of Beam
    material = steel,   # Material
    section = section   # Crossection
)

# __ Discretisation and Meshing _______
n = 20
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
Load1 = lds.PointLoad(nodemap[loadpoint], Fy=-100)

U, _ = sol.static(loading=[Load1], boundaryconditions=[simplesupport1, simplesupport2], K=K)
print(f"The deflection of simply supported beam at middle using FEM = {U[Load1.node.ID*3 + 1]}")


umax = (Load1.Fy * (beam1.length() ** 3)) / (48 * steel.E * section.Izz)
print(f"The deflection of simply supported beam at middle with analytical solutions = {umax}")

err = (np.abs((U[Load1.node.ID*3 + 1] - umax) / umax)) * 100
print(f"Percentage error = {err}")

pstprcs.EBBTpostprocess(U, elements)
