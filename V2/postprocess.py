import numpy as np

import V2.entities as ent

def EBBTpostprocess(U, elements, K, F):
    def kinematics(element: ent.element1D, u_l, x = None):
        def dNa_dx(x, L):
            dN1_dx = -1/L
            dN2_dx =  1/L
            return np.array([dN1_dx, dN2_dx])
            
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
        
        dl = element.length()
        h = element.parent[0].section.thickness
        if x is None: x = np.linspace(0, dl, 6)
        
        u = u_l[[0, 3]]
        v = u_l[[1, 2, 4, 5]]

        v_xx = (d2Nb_dx2(x, dl).T @ v)
        v_xxx = (d3Nb_dx3(x, dl).T @ v)
        u_x = (np.tile(dNa_dx(x, dl).reshape(-1, 1), (1, 6)).T @ u)

        return {
            "eps0" : u_x,
            "kappa" : v_xx,
            "d3 v/dx3" : v_xxx,
        }

    def element_disp(U, element : ent.element1D):
        return element.T() @ U[element.dof()]
        
    def strain(element: ent.element1D, u_l):
        diffs = kinematics(element, u_l)

        u_x = diffs["eps0"]
        v_xx  = diffs["kappa"]

        h = element.parent[0].section.thickness

        y = np.linspace(-h/2, h/2, 5)
        eps = u_x + y[:, None]*v_xx
        
        return eps

    def stress(element: ent.element1D, strain):
        return (element.parent[0].material.E)*strain
    
    def reactions(K, U, F):
        return (K @ U - F)

    def internalforces(element: ent.element1D, u_l):
        diffs = kinematics(element, u_l)
        
        E = element.parent[0].material.E
        I = element.parent[0].section.Izz
        A = element.parent[0].section.A

        u_x = diffs["eps0"]
        v_xx = diffs["kappa"]
        v_xxx = diffs["d3 v/dx3"]

        N = E * A * u_x
        M = E * I * v_xx
        V = E * I * v_xxx

        return {
            "axial" : N,
            "moment" : M,
            "shear" : V
        }
    
    Strains = []
    Stresses = []
    Reactions = reactions(K, U, F)
    Internalforces = []
    for element in elements:
        ul = element_disp(U, element)
        eps = strain(element, ul)
        Strains.append(eps)
        Stresses.append(stress(element, eps))
        Internalforces.append(internalforces(element, ul))

    return{
        "Strain" : Strains,
        "Stress" : Stresses,
        "Reactions" : Reactions,
        "Internal Forces" : Internalforces,
    }
    
