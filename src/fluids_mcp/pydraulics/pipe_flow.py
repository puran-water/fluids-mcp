import fluids

from math import pi
#Constants

G = 9.81 #gravity
nu_20c = 0.0000010533 #kinematic viscosity at 20 C

class Pipe:
    def __init__(self, Q, D, e= None, L= None, rho = None, mu = None, C = None, method= 'Darcy-Weisbach'):
        self.Q = Q
        self.D = D
        self.a = self.D**2 * pi / 4
        self.e = e
        self.eD = self.e / self.D if (e and D) != None else None
        self.v = Q / (pi * D**2 /4)
        self.nu = rho / mu if (rho and mu) != None else nu_20c
        self.Re = self.v * self.D / self.nu
        self.L = L
        self.method = method
        self.C = C
        self.h = self.calc_head_loss()
        self.hf = self.h * self.L if L is not None else None
        self.h = round(self.h, 4) if self.h is not None else None
        self.hf = round(self.hf, 4) if self.hf is not None else None


    def calc_head_loss(self):
        match self.method:
            case 'Darcy-Weisbach':
                def calc_friction(self):
                   return fluids.friction_factor(self.Re, self.eD)
                if(self.Re < 2000):
                    self.flow_type = 'Laminar flow'
                    self.fr = 64 / self.Re
                elif (self.Re < 4000):
                    self.flow_type = 'Transitional flow'
                    self.fr = calc_friction(self)
                elif (self.Re >= 4000):
                    self.flow_type = 'Turbulent flow'
                    self.fr = calc_friction(self)
                if self.L is not None:
                    h = self.fr * self.v**2 / (2 * self.D * G)
                    self.fr = round(self.fr, 4)
                    return h
                else:
                    self.fr = round(self.fr, 4) if hasattr(self, 'fr') else None
                    return None
            case 'Hazen-Williams':
                if self.L is not None:
                    h = 10.674 * self.Q**1.852 / (self.D**4.87 * self.C**1.852)
                    return h
                else:
                    return None
