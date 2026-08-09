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

section1 = sec.Circle(radius = 0.01)
steel = mat.Material()

beam1 = defin.CurvedBeam_element(start = (-1, 0), end = (0, 1), radius = 1, material=steel, section=section1)
print(beam1.center())
print(np.rad2deg(beam1.start_angle()), np.rad2deg(beam1.end_angle()))
beam2 = defin.StraightBeam_element(start= (0, 1), end= (0, -2), material=steel, section=section1)
beam3 = defin.StraightBeam_element(start= (-1, 0), end= (0, 0), material=steel, section=section1)
beam4 = defin.StraightBeam_element(start= (-1, 0), end= (-1, -2), material=steel, section=section1)
beam5 = defin.StraightBeam_element(start= (-1, -2), end= (0, -2), material=steel, section=section1)


nodes, elements, nodesmap = mesh.OneDMeshing(np.array([(beam1, 300), (beam2, 300), (beam3, 100), (beam4, 200), (beam5, 200)]))

print(len(nodes), len(elements))

plt.figure()
for element in elements:
    x = [element.node1.coords[0], element.node2.coords[0]]
    y = [element.node1.coords[1], element.node2.coords[1]]
    plt.plot(x, y, color = (0, 1, 1))
plt.axis("equal")
plt.grid()
plt.title("System underformed")
plt.savefig("arctesting/undeformedstate.png")
plt.show()

K, M = ass.assemble(nodes, elements)

loadpoint = (0, 1)
load = lds.PointLoad(nodesmap[(0, 1)], Fy= -1000000)
load2 = lds.PointLoad(nodesmap[(0, 0)], Fx=-1000)

pinsupport1 = lds.PinSupport(nodesmap[(-1, -2)])
pinsupport2 = lds.PinSupport(nodesmap[(0, -2)])

staticresults = sol.static(loading=[load, load2], boundaryconditions=[pinsupport1, pinsupport2], K=K)
deflections = staticresults["Deformation"]

plt.figure()
for element in elements:
    dofs = element.dof()
    u = deflections[dofs[[0, 3]]]
    v = deflections[dofs[[1, 4]]]
    x = [element.node1.coords[0], element.node2.coords[0]] + u
    y = [element.node1.coords[1], element.node2.coords[1]] + v
    plt.plot(x, y, color = (1, 0, 1))
plt.axis("equal")
plt.grid()
plt.title("System deformed under a load of -10000 vertical load \nat (0, 1)")
plt.savefig("arctesting/deformedstate1.png")
plt.show()