# __ DEFINING NODES ________________
class node:
    def __init__(self, ID, coords):
        self.ID = ID            # Unique ID for each Node
        self.coords = coords    # Coordinates of the Node
        
# __ ELEMENT ENTITIES ______________
class element1D:
    def __init__(self, elementID, node1, node2, parent, extranodes = None):
        self.ID = elementID         # Unique ID for each FE element
        self.node1 = node1          # Start Node of the element
        self.node2 = node2          # End Node of the element 
        self.parent = parent        # The Original beam object it came from 
        self.nextnodes = extranodes # For p-refinement


