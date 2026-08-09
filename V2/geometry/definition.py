import numpy as np

# __ GEOMETRY DEFINITION ______________
class StraightBeam_element:
    def __init__(self, start, end, material, section):
        self.start = start
        self.end = end
        self.material = material
        self.section = section
        self.type = "Straight"

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


class CurvedBeam_element:
    def __init__(self, start, end, radius, material, section, major_arc = False, ):
        # Clockwise direction is used as convention for the start and end points of the curved beam.
        # The start and end points are defined such that the beam follows a clockwise path.
        self.start = start              # Start Point of curved Beam
        self.end = end                  # End Point of curved Beam
        self.radius = radius            # Radius of curved Beam
        self.material = material        # Material of curved Beam
        self.section = section          # Section of curved Beam
        self.major_arc = major_arc      # True if the curved beam is a major arc, False if it is a minor arc
        self.type = "Curved"

    def length(self):
            S = self.start
            E = self.end
            return np.sqrt((E[0] - S[0]) ** 2 + (E[1] - S[1]) ** 2)

    def center(self):
        S = self.start
        E = self.end
        R = self.radius
        L = self.length()
        if L > 2 * R:
            raise ValueError("The length of the chord is greater than the diameter of the circle.")
        else:
            if((self.major_arc)):
                x_c1 = (S[0] + E[0]) / 2 + np.sqrt(R ** 2 - (L / 2) ** 2) * (S[1] - E[1]) / L
                y_c1 = (S[1] + E[1]) / 2 + np.sqrt(R ** 2 - (L / 2) ** 2) * (E[0] - S[0]) / L
                            
                return (x_c1, y_c1)
            
            elif(not self.major_arc):
                x_c2 = (S[0] + E[0]) / 2 - np.sqrt(R ** 2 - (L / 2) ** 2) * (S[1] - E[1]) / L
                y_c2 = (S[1] + E[1]) / 2 - np.sqrt(R ** 2 - (L / 2) ** 2) * (E[0] - S[0]) / L
                
                return (x_c2, y_c2)

    
    def start_angle(self):
        S = self.start
        C = self.center()
        theta = np.atan2((S[1] - C[1]), (S[0] - C[0]))
        theta = theta + 2*np.pi if theta < 0 else theta

        return theta

    def end_angle(self):
        E = self.end
        C = self.center()
        theta = np.atan2((E[1] - C[1]), (E[0] - C[0]))
        theta = theta + 2 * np.pi if theta < 0 else theta

        return theta

    def arc_length(self):
        start_angle = self.start_angle()
        end_angle = self.end_angle()

        angle_diff = np.abs(end_angle - start_angle)
        if self.major_arc:
            angle_diff = 2 * np.pi - angle_diff

        return self.radius * angle_diff


class Spring_element:
    def __init__(self, first, second, springconstant, mass1 = 0, mass2 = 0):
        E = second
        S = first
        self.first = S
        self.second = E
        self.springconstant = springconstant
        self.naturallength = np.sqrt((E[0] - S[0]) ** 2 + (E[1] - S[1]) ** 2)
        self.mass1 = mass1
        self.mass2 = mass2
