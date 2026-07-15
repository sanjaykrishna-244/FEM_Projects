import numpy as np
import V2.entities as ent

# __ BOUNDARY CONDITOINS ______________
class FixedSupport:
    def __init__(self, node:ent.node):
        self.node = node

    def fixeddof(self):
        return [self.node.ID * 3, self.node.ID * 3 + 1, self.node.ID * 3 + 2]

class PinSupport:
    def __init__(self, node:ent.node):
        self.node = node
    
    def fixeddof(self):
        return [self.node.ID * 3, self.node.ID * 3 + 1]

class RollerSupport:
    def __init__(self, node:ent.node, direction: str):
        self.node = node
        self.direction = direction      # the direction constrained to move
    
    def fixeddof(self):
        if self.direction == "x":
            return [self.node.ID * 3]
        elif self.direction == "y":
            return [self.node.ID * 3 + 1]
        

# __ LOADS ____________________________
class PointLoad:
    def __init__(self, node:ent.node, Fx:float = 0.0, Fy:float = 0.0, M:float = 0.0):
        self.node = node
        self.Fx = Fx
        self.Fy = Fy
        self.M = M
        self.type = "pointload"

