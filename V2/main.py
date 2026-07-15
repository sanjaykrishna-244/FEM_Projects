import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import V2.geometry.crosssection as sec
import V2.geometry.definition as define
import V2.material as mat

import V2.mesh as mesh
import V2.asssemble as ass

import V2.loads as lds

import V2.solver as sol

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
n = 40
frame = np.array([(beam1, n), (beam2, n), (beam3, n), (beam4, n), (beam5, n), (beam6, n), (beam7, n)])
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
fixedsupportpoint1 = (0, 0)
fixedsupport1 = lds.FixedSupport(nodesmap[fixedsupportpoint1])

rollersupportpoint = (1, 0)
rollersupport = lds.RollerSupport(nodesmap[rollersupportpoint], "y")

pinsupportpoint = (2, 0)
pinsupport = lds.PinSupport(nodesmap[pinsupportpoint])

print(rollersupport.fixeddof())

loadpoint = (1, 1)
Load1 = lds.PointLoad(nodesmap[loadpoint], Fy = -2000)

N = K.shape[0]

U = sol.static(loads = [Load1], boundaryconditions=[fixedsupport1, pinsupport], K = K)
print(U[Load1.node.ID*3 + 1])
Ux = (U[np.arange(0, N, 3)] + XY[:, 0]).reshape(-1, 1)
Uy = (U[np.arange(1, N, 3)] + XY[:, 1]).reshape(-1, 1)
U_ = np.hstack((Ux, Uy, IDs))
print(U_.shape)

#####
plt.figure(figsize=(8,4))

plt.scatter(U_[:,0], U_[:,1], color='k')

for (x, y, id) in U_:
    plt.text(x, y, str(int(id)),
             fontsize=8,
             color='red')

plt.axis("equal")
plt.grid(True)
plt.show()


f, modes = sol.modal(boundaryconditions=[fixedsupport1, pinsupport], K = K, M = M)
print(f[0])
U = modes[:, 4]
Ux = (0.5*U[np.arange(0, N, 3)] + XY[:, 0]).reshape(-1, 1)
Uy = (0.5*U[np.arange(1, N, 3)] + XY[:, 1]).reshape(-1, 1)
U_ = np.hstack((Ux, Uy, IDs))
print(U_.shape)

#####
plt.figure(figsize=(8,4))

plt.scatter(U_[:,0], U_[:,1], color='k')

for (x, y, id) in U_:
    plt.text(x, y, str(int(id)),
             fontsize=8,
             color='red')

plt.axis("equal")
plt.grid(True)
plt.show()


'''
dfk = pd.DataFrame(K)
dfm = pd.DataFrame(M)

with pd.ExcelWriter("sample.xlsx", engine = 'openpyxl') as writer:
    dfk.to_excel(writer, sheet_name="stiffness", index=False)
    dfm.to_excel(writer, sheet_name="Mass", index=False)
'''