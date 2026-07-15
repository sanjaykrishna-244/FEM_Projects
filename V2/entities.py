import numpy as np

# __ DEFINING NODES ________________
class node:
    def __init__(self, ID:int, coords:tuple):
        self.ID = ID            # Unique ID for each Node
        self.coords = coords    # Coordinates of the Node
        
# __ ELEMENT ENTITIES ______________
class element1D:
    def __init__(self, elementID:int, node1:node, node2:node, parent, extranodes = None):
        self.ID = elementID         # Unique ID for each FE element
        self.node1 = node1          # Node1 of the element
        self.node2 = node2          # Node2 of the element 
        self.parent = parent        # The Original beam object it came from 
        self.nextnodes = extranodes # For p-refinement

    def length(self):
        node1 = self.node1.coords
        node2 = self.node2.coords
        return np.round(np.sqrt((node2[0] - node1[0])**2 + (node2[1] - node1[1])**2), 8)
    
    def angle(self):
        node1 = self.node1.coords
        node2 = self.node2.coords
        return np.round(np.atan2((node2[1] - node1[1]), (node2[0] - node1[0])), 8)
    
    def dof(self):
        id1 = self.node1.ID
        id2 = self.node2.ID
        dofs = np.array([3*id1, 3*id1 + 1, 3*id1 + 2, 3*id2, 3*id2 + 1, 3*id2 + 2])
        return dofs
    
    def T(self):
        c = np.round(np.cos(self.angle()), 5)
        s = np.round(np.sin(self.angle()), 5)

        T = np.array([[ c, s, 0, 0, 0, 0],
                      [-s, c, 0, 0, 0, 0],
                      [ 0, 0, 1, 0, 0, 0],
                      [ 0, 0, 0, c, s, 0],
                      [ 0, 0, 0,-s, c, 0],
                      [ 0, 0, 0, 0, 0, 1]])

        return T
    
    def Ke(self):
        dl = self.length()
        T = self.T()
        E = self.parent[0].material.E
        A = self.parent[0].section.A
        I = self.parent[0].section.Izz
        a = 12; b = 6 * dl; c = 2 * dl ** 2
        k = np.array([
                    [1, 0, 0, -1, 0, 0],
                    [0, a, b, 0, -a, b], 
                    [0, b, 2*c, 0, -b, c], 
                    [-1, 0, 0, 1, 0, 0], 
                    [0, -a, -b, 0, a, -b], 
                    [0, b, c, 0, -b, 2*c]
                    ])
        axial = [0, 3]
        bendg = np.setdiff1d(np.arange(0, 6), axial)
        k[np.ix_(axial, axial)] *= E*A/dl
        k[np.ix_(bendg, bendg)] *= E*I/(dl**3)
        k = np.transpose(T) @ k @ T

        return k
    
    def Me(self):
        dl = self.length()
        T = self.T()
        rho = self.parent[0].material.rho
        A = self.parent[0].section.A
        m = np.array([
            [140, 0, 0, 70, 0, 0],
            [0, 156, 22*dl, 0, 54, -13*dl],
            [0, 22*dl, 4*dl**2, 0, 13*dl, -3*dl**2],
            [70, 0, 0, 140, 0, 0],
            [0, 54, 13*dl, 0, 156, -22*dl],
            [0, -13*dl, -3*dl**2, 0, -22*dl, 4*dl**2]]
            )
        m *= rho * A * dl / 420
        m = T.T @ m @ T

        return m