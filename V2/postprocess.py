import numpy as np

import V2.entities as ent

def EBBTpostprocess(U, elements):
    def element_disp(U, element : ent.element1D):
        return element.T() @ U[element.dof()]
        
    def strain(element: ent.element1D, u_l):
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
        
        dl = element.length()
        N = 5
        h = element.parent[0].section.thickness
        x = np.arange(0, N+1) * (dl / N)
        y = np.linspace(-h/2, h/2, 5)
        u = u_l[[0, 3]]
        v = u_l[[1, 2, 4, 5]]

        v_xx = (d2Nb_dx2(x, dl).T @ v)
        u_x = (np.tile(dNa_dx(x, dl).reshape(-1, 1), (1, 6)).T @ u)

        eps = u_x + y[:, None]*v_xx

        return eps

    def stress(element: ent.element1D, strain):
        return (element.parent[0].material.E)*strain
    
    def reactions(K, U, F):
        return (K @ U - F)
    
     
    element1 = elements[0]
    eps1 = strain(element1, element_disp(U, element1))
    print(stress(element1, eps1))
    

