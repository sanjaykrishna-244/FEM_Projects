import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import V2.geometry.crosssection as sec
import V2.geometry.definition as define
import V2.material as mat

import V2.mesh as mesh

import V2.asssemble as ass

steel = mat.Material()
rect = sec.Rectangle()
beam1 = define.Beam_element((0, 0), (1, 0), steel, rect)
beam2 = define.Beam_element((1, 0), (2, 0), steel, rect)
beam3 = define.Beam_element((0, 1), (1, 1), steel, rect)
beam4 = define.Beam_element((1, 1), (2, 1), steel, rect)
beam5 = define.Beam_element((0, 0), (0, 1), steel, rect)
beam6 = define.Beam_element((1, 0), (1, 1), steel, rect)
beam7 = define.Beam_element((2, 0), (2, 1), steel, rect)
#print(beam3.sin())

frame = np.array([(beam1, 20), (beam2, 20), (beam3, 20), (beam4, 20), (beam5, 20), (beam6, 20), (beam7, 20)])
nodes, elements, nodesmap = mesh.OneDMeshing(frame)

XY = np.array([list(node.coords) for node in nodes])
IDs = np.array([node.ID for node in nodes]).reshape(-1, 1)
XY = np.hstack((XY, IDs))
#print(XY[7])
#print(nodes[0].ID, nodes[0].coords)
#print(elements[2].ID, elements[2].node1.ID, elements[2].node2.ID)
#print(np.rad2deg(elements[90].angle()))
#print(elements[100].Ke())

K, M = ass.assemble(nodes, elements)
#print(K.shape)
#req = elements[100].dof()
#print(elements[100].node1.ID)
#print(elements[100].node2.ID)

#print(K[np.ix_(req, req)])
#####
plt.figure(figsize=(8,4))

plt.scatter(XY[:,0], XY[:,1], color='k')

for (x, y, id) in XY:
    plt.text(x, y, str(int(id)),
             fontsize=8,
             color='red')

plt.axis("equal")
plt.grid(True)
plt.show()

dfk = pd.DataFrame(K)
dfm = pd.DataFrame(M)

with pd.ExcelWriter("sample.xlsx", engine = 'openpyxl') as writer:
    dfk.to_excel(writer, sheet_name="stiffness", index=False)
    dfm.to_excel(writer, sheet_name="Mass", index=False)
