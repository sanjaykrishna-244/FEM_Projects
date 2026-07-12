from numpy import pi as pi

# __ CROSS-SECTION DEFINTIONS _________
# __ Rectangular cross-section __
class Rectangle:
    def __init__(self, width = 0.01, thickness = 0.01):
        self.width = width
        self.thickness = thickness
        self.A = thickness * width
        self.Izz = (width) * (thickness ** 3.0) / 12.0
        self.Iyy = (width ** 3.0) * (thickness) / 12.0

# __ Circular cross-section __
class Circle:
    def __init__(self, radius = 0.01):
        self.radius = radius
        self.A = pi * (radius ** 2.0)
        self.I = pi * (radius ** 4.0) / 4.0