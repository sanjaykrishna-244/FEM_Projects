import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import openpyxl as op
import utlities as ut
from shipdeckframe.load1 import u
from shipdeckframe.doubleframe import details, E, t, w


def Na(x, L):
    N1 = 1 - (x / L)
    N2 = (x / L)
    return np.array([N1, N2])

def dNa_dx(x, L):
    dN1_dx = -1/L
    dN2_dx =  1/L
    return np.array([dN1_dx, dN2_dx])

def Nb(x, L):
    H1 = 1 - 3*(x/L)**2 + 2*(x/L)**3
    H2 = x - 2*(x**2)/L + (x**3)/(L**2)
    H3 = 3*(x/L)**2 - 2*(x/L)**3
    H4 = - (x**2)/L + (x**3)/(L**2)
    return np.array([H1, H2, H3, H4])

def dNb_dx(x, L):
    dH1_dx = - 6*(x)/(L**2) + 6*(x**2)/(L**3)
    dH2_dx = 1 - 4*(x)/L + 3*(x**2)/(L**2)
    dH3_dx = 6*(x)/(L**2) - 6*(x**2)/(L**3)
    dH4_dx = -2*(x)/L + 3*(x**2)/(L**2)
    return np.array([dH1_dx, dH2_dx, dH3_dx, dH4_dx])

def d2Nb_dx2(x, L):
    d2H1_dx2 = - 6/(L**2) + 12*(x)/(L**3)
    d2H2_dx2 = - 4/L + 6*x/(L**2)
    d2H3_dx2 = 6/(L**2) - 12*x/(L**3)
    d2H4_dx2 = -2/L + 6*x/(L**2)
    return np.array([d2H1_dx2, d2H2_dx2, d2H3_dx2, d2H4_dx2])

def d3Nb_dx3(x, L):
    d3H1_dx3 = 12/(L**3)
    d3H2_dx3 = 6/(L**2)
    d3H3_dx3 = - 12/(L**3)
    d3H4_dx3 = 6/(L**2)
    return np.array([d3H1_dx3, d3H2_dx3, d3H3_dx3, d3H4_dx3])

def stressextract(dofs, u, details, id):
    u_l = u[dofs]
    T = ut.transformation(details[id][1])
    u_l = T @ u_l
    x = np.arange(0, 6) * details[id][2] / 5
    uvals = u_l[[0, 3]]
    vvals = u_l[[1, 2, 4, 5]]
    U = np.zeros_like(x)
    V = np.zeros_like(x)
    for i in range(len(x)):
        N = dNa_dx(x[i], details[id][2])
        B = d2Nb_dx2(x[i], details[id][2])
        u_x  = N @ uvals
        v_xx = B @ vvals
        U[i] = u_x
        V[i] = v_xx
    y = np.arange(-2, 3) * t / 4
    # print(U, V)
    U = np.repeat(U[np.newaxis, :], repeats=y.shape[0], axis = 0)
    y = y.reshape(len(y), 1)
    V = V.reshape(1, len(V))
    strain = U - (y @ V)
    stress = E * strain
    return stress

elems = np.array([[1, 2, 0, 0, 0]])
for i in details:
    node1s = details[i][0][:-1]
    node2s = details[i][0][1:]
    ids = i*np.ones(len(node1s), dtype = int)
    x1s = details[i][3][:-1]
    y1s = details[i][4][:-1] 
    elem = np.array([node1s, node2s, ids, x1s, y1s]).T
    elems = np.concatenate([elems, elem], axis = 0)
elems = elems[1:, :]
#print(nodes)

framestress = []
for node in elems:
    N1 = int(node[0])
    N2 = int(node[1])
    dofs = np.array([3*N1, 3*N1+1, 3*N1+2, 3*N2, 3*N2+1, 3*N2+2])
    stress = stressextract(dofs, u, details, node[2])
    framestress.append(stress)

#print(framestress[85])

x = np.arange(0,6) * 0.04
y = np.arange(-2,3) * t / 4
grid = np.array(np.meshgrid(x, y))
grid = np.moveaxis(grid, 0, -1)

finalcoords = np.array([])
for i in range(len(elems)):
    T = ut.coordinate_transfor(details[elems[i][2]][1])
    coords = grid @ T + np.array([elems[i][3], elems[i][4]])
    coords = np.round(coords, 3)
    if i == 0:
        finalcoords = np.array([coords])
    else:
        finalcoords = np.append(finalcoords, [coords], axis = 0)

fig, ax = plt.subplots()
for e in range(len(elems)):
    x = finalcoords[e, :, :, 0]
    y = finalcoords[e, :, :, 1]
    stress = framestress[e]
    pcm = ax.pcolormesh(
        x,
        y,
        stress,
        shading='gouraud',
        cmap='RdBu_r'
    )
fig.colorbar(pcm, ax=ax, label="Stress (Pa)")
plt.axis("equal")
plt.show()