# __ ISOTROPIC MATERIAL PROPERTIES ____
class Material:
    def __init__(self, E = 2e11, nu =0.3, rho = 7850, sigma_y = 355e6, sigma_u = 475e6):
        self.E = E
        self.nu = nu
        self.rho = rho
        self.sigma_y = sigma_y
        self.sigma_u = sigma_u

    def G(self):
        E = self.E
        n = self.nu
        return E/(2*(1+n))

    def strain(self, stress):
        if stress < self.sigma_y:
            return stress / self.E
        return -1  

# __ PLASTICITY _______________________
