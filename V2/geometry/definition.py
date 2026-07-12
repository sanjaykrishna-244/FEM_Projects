import numpy as np

# __ GEOMETRY DEFINITION ______________
class Beam_element:
    def __init__(self, start, end, material, section):
        self.start = start
        self.end = end
        self.material = material
        self.section = section

    def length(self):
        S = self.start
        E = self.end
        return np.sqrt((E[0] - S[0]) ** 2 + (E[1] - S[1]) ** 2)
    
    def angle(self):
        S = self.start
        E = self.end
        return np.atan2((E[1] - S[1]), (E[0] - S[0]))
    
    def cos(self):
        return (self.end[0] - self.start[0]) / self.length()
    
    def sin(self):
        return (self.end[1] - self.start[1]) / self.length()


