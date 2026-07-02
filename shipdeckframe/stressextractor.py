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

def stressextract(dofs, u, details):
    u_l = u[dofs[1]]
    T = ut.transformation(details[1][1])
    u_l = T @ u_l
    x = np.arange(0, 11) * 0.02
    uvals = u_l[[0, 3]]
    vvals = u_l[[1, 2, 4, 5]]
    U = np.zeros_like(x)
    V = np.zeros_like(x)
    for i in range(len(x)):
        N = dNa_dx(x[i], details[1][2])
        B = d2Nb_dx2(x[i], details[1][2])
        u_x  = N @ uvals
        v_xx = B @ vvals
        U[i] = u_x
        V[i] = v_xx
    y = np.arange(-2, 3) * t / 4
    U = np.repeat(U[np.newaxis, :], repeats=y.shape[0], axis = 0)
    y = y.reshape(len(y), 1)
    V = V.reshape(1, len(V))
    stress = U + (y @ V)
    
    return stress

element_nodes = [0, 1]
dof = {1 : np.array([3*element_nodes[0], 3*element_nodes[0]+1, 3*element_nodes[0]+2, 3*element_nodes[1], 3*element_nodes[1]+1, 3*element_nodes[1]+2])}

ele1stress = stressextract(dof, u, details)
ele1stress = np.repeat(ele1stress[np.newaxis, :, :], repeats=9, axis = 0)
print(ele1stress.shape)