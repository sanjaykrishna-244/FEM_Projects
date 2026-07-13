import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import V2.geometry.crosssection as sec
import V2.geometry.definition as define
import V2.material as mat

import V2.mesh as mesh

steel = mat.Material()
rect = sec.Rectangle()
beam1 = define.Beam_element((0, 0), (1, 0), steel, rect)
beam2 = define.Beam_element((1, 0), (1, 1), steel, rect)
beam3 = define.Beam_element((1, 1), (0, 0), steel, rect)

print(beam3.sin())

frame = np.array([(beam1, 20), (beam2, 20), (beam3, 30)])

nodes, elements = mesh.OneDMeshing(frame)

XY = np.array([list(node.coords) for node in nodes])
IDs = np.array([node.ID for node in nodes]).reshape(-1, 1)
XY = np.hstack((XY, IDs))
print(XY[7])
print(nodes[0].ID, nodes[0].coords)
print(elements[2].ID, elements[2].node1.ID, elements[2].node2.ID)


plt.figure(figsize=(8,4))

plt.scatter(XY[:,0], XY[:,1], color='k')

for (x, y, id) in XY:
    plt.text(x, y, str(int(id)),
             fontsize=8,
             color='red')

plt.axis("equal")
plt.grid(True)
plt.show()
