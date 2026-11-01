# Imports
from sympy import Symbol, nsolve, sqrt, acos, sin, lambdify, re
import numpy as np 
import sys


# Variables definition

'''
https://www.eng.auburn.edu/~xzf0001/Handbook/Channels.html

'n': "Manning's Coef.",
'So': 'Channel Slope [m/m]',
'Q': 'Flow Rate [m3/s]',
'y': 'Depth [m]',
'b': 'Bottom Width [m]',
'z': 'Side Slope',
'D': 'Diameter [m]',
'a': 'Area [m2]',
'rh': 'Hydraulic Radius',
'dh': 'Hydraulic Depth',
'tw': 'Top Width [m]',
'p': 'Wetted Perimeter [m]',
'f': 'Froude Number',
'v': 'Velocity [m/s]',
'flow_status': 'Flow Status',
'zc': 'Section Factor',
'yc': 'WIP',

'''

# Constants

G = 9.81 #acceleration due gravity

# Classes

class Channel:
    def __init__(self, n, So, Q):
        self.n = n
        self.So = So
        self.Q = Q
        self.y = None
        self.a = None
        self.p = None
        self.rh = None
        self.tw = None
        self.dh = None
        self.zc = None
        self.v = None
        self.f = None
        self.flow_status = None
        self.yn = Symbol('yn')
        self.ac = None
        self.rc = None
        self.pc = None
        self.twc = None
        self.yc = None
        self.Sc = None
        
    
    def calc_properties(self):
        if self.Q is not None and type(self.y) is not Symbol:
            # Only calculate if prerequisites are valid
            if self.a is not None and self.a > 0:
                self.v = self.Q / self.a
                if self.dh is not None and self.dh > 0:
                    try:
                        self.f = self.v / sqrt(G * self.dh)
                        froude_value = float(str(self.f))
                        if froude_value > 1.0:
                            self.flow_status = 'Supercritical'
                        elif abs(froude_value - 1.0) < 1e-6:
                            self.flow_status = 'Critical'
                        else:
                            self.flow_status = 'Subcritical'
                    except (ValueError, TypeError):
                        self.f = None
                        self.flow_status = None
            
            # Calculate the hydraulic radius from critical values if available
            if self.ac is not None and self.pc is not None and self.pc > 0:
                self.rc = self.ac / self.pc
            else:
                self.rc = None
        else:
            pass

    def calc_flow(self):
        #Manning
        self.calc_properties()
        self.Q = (self.a * self.rh**(2/3) * self.So**0.5) / self.n
        self.calc_properties()
        return self.Q

    def calc_yn(self):
        if(type(self.y) is Symbol):
            try:
                self.calc_properties()
                sol = (self.a * self.rh**(2/3) * self.So**0.5) / self.n
                self.y = re(nsolve(sol - self.Q, self.y, 1))
                self.calc_properties()
            except ValueError:
                return 'Incorrect value selected'
        return self.y

    def get_parameters(self):
        return self.__dict__
    
    def get_energy(self):
        y = np.arange(0.1,4, 0.01) 
        # Crear una función evaluable a partir de la expresión self.ac
        specific_energy = lambdify(self.yn, self.yn + (self.Q**2 / (2 * G * self.ac**2)))
        # Evaluar la función en los puntos de y
        E_values = specific_energy(y)
        return y, E_values

    def get_critical_parameters(self):
        froude_critical = (self.Q**2 * self.twc) - (self.ac**3 * G)
        x0 = 1 if (hasattr(self, 'D') and self.D > 1.1) else 0.4
        self.yn = re(nsolve(froude_critical, self.yn, x0))
        self.yc = self.yn
        self.calc_properties()
        self.Sc = ((self.Q**2 * self.n**2) / (self.ac**2 * self.rc**(4/3)))


class RectangularChannel(Channel):
    def __init__(self, n, So, b, Q = None, y = Symbol('y')):
        super().__init__(n, So, Q)
        self.channel_type = 'Rectangular'
        self.b = b
        self.y = y if y != None else Symbol('y')
        if(type(self.y) is Symbol):
            self.calc_yn()
        if(self.Q is None):
            self.calc_flow()
        self.get_critical_parameters()
        

    def calc_properties(self):
        try:
            # Handle numeric y calculations
            if self.y is not None and type(self.y) is not Symbol:
                try:
                    self.a = float(self.b) * float(self.y)
                    self.p = float(self.b) + (float(self.y) * 2)
                    if self.p > 0:
                        self.rh = self.a / self.p
                    else:
                        self.rh = None
                    self.tw = float(self.b)
                    self.dh = float(self.y)
                    self.zc = float(self.b) * (float(self.y)**1.5)
                except (ValueError, TypeError):
                    self.a = self.p = self.rh = self.tw = self.dh = self.zc = None
            
            # Handle symbolic yn calculations
            if self.yn is not None:
                try:
                    self.ac = self.b * self.yn
                    self.pc = self.b + (self.yn * 2)
                    self.twc = float(self.b)
                except (ValueError, TypeError):
                    self.ac = self.pc = self.twc = None
            
            # Call parent method only if essential values are set
            if hasattr(self, 'a') and self.a is not None:
                super().calc_properties()
        except Exception as e:
            # Handle any unexpected errors
            self.a = self.p = self.rh = self.tw = self.dh = self.zc = None

    
    
class TrapezoidalChannel(Channel):
    def __init__(self, n, So, z, b, Q = None, y = Symbol('y')):
        super().__init__(n, So, Q)
        self.channel_type = 'Trapezoidal'
        self.b = b
        self.z = z
        self.y = y if y != None else Symbol('y')
        if(type(self.y) is Symbol):
            self.calc_yn()
        if(self.Q is None):
            self.calc_flow()
        self.get_critical_parameters()
        

    def calc_properties(self):
        try:
            # Handle numeric y calculations
            if self.y is not None and type(self.y) is not Symbol:
                try:
                    b = float(self.b)
                    y = float(self.y)
                    z = float(self.z)
                    
                    self.a = (b + (y * z)) * y
                    self.p = b + (2 * y * (1 + z**2)**(1/2))
                    if self.p > 0:
                        self.rh = self.a / self.p
                    else:
                        self.rh = None
                        
                    self.tw = b + (2 * y * z)
                    
                    tw_denom = (b + (2 * z * y))
                    if tw_denom > 0:
                        self.dh = ((b + (z * y)) * y) / tw_denom
                    else:
                        self.dh = None
                        
                    zc_denom = (b + (2 * y * b))
                    if zc_denom > 0:
                        self.zc = ((b + (z * y)) * y)**1.5 / zc_denom**0.5
                    else:
                        self.zc = None
                except (ValueError, TypeError):
                    self.a = self.p = self.rh = self.tw = self.dh = self.zc = None
            
            # Handle symbolic yn calculations
            if self.yn is not None:
                try:
                    self.ac = (self.b + (self.yn * self.z)) * self.yn
                    self.twc = self.b + (2 * self.yn * self.z)
                    self.pc = self.b + (2 * self.yn * (1 + self.z**2)**(1/2))
                except (ValueError, TypeError):
                    self.ac = self.twc = self.pc = None
            
            # Call parent method only if essential values are set
            if hasattr(self, 'a') and self.a is not None:
                super().calc_properties()
        except Exception as e:
            # Handle any unexpected errors
            self.a = self.p = self.rh = self.tw = self.dh = self.zc = None


class TriangularChannel(Channel):
    def __init__(self, n, So, z, Q = None, y = Symbol('y')):
        super().__init__(n, So, Q)
        self.channel_type = 'Triangular'
        self.z = z
        self.y = y if y != None else Symbol('y')
        if(type(self.y) is Symbol):
            self.calc_yn()
        if(self.Q is None):
            self.calc_flow()
        self.get_critical_parameters()
        

    def calc_properties(self):
        try:
            # Handle numeric y calculations
            if self.y is not None and type(self.y) is not Symbol:
                try:
                    y = float(self.y)
                    z = float(self.z)
                    
                    self.a = z * y**2
                    self.p = (y * 2) * (1 + z**2)**0.5
                    if self.p > 0:
                        self.rh = self.a / self.p
                    else:
                        self.rh = None
                        
                    self.tw = 2 * z * y
                    self.dh = 0.5 * y
                    self.zc = ((2**0.5)/2) * z * y**2.5
                except (ValueError, TypeError):
                    self.a = self.p = self.rh = self.tw = self.dh = self.zc = None
            
            # Handle symbolic yn calculations
            if self.yn is not None:
                try:
                    self.an = self.z * self.yn**2
                    self.ac = self.z * self.yn**2
                    self.twc = 2 * self.z * self.yn
                    self.pc = (self.yn * 2) * (1 + self.z**2)**0.5
                except (ValueError, TypeError):
                    self.an = self.ac = self.twc = self.pc = None
            
            # Call parent method only if essential values are set
            if hasattr(self, 'a') and self.a is not None:
                super().calc_properties()
        except Exception as e:
            # Handle any unexpected errors
            self.a = self.p = self.rh = self.tw = self.dh = self.zc = None

    def caudal_manning(self):
        self.calc_properties()
        self.Q = (self.a * self.rh**(2/3) * self.So**0.5) / self.n
        self.calc_properties()    
    

class CircularChannel(Channel):
    def __init__(self, n, So, D, Q = None, y = Symbol('y')):
        super().__init__(n, So, Q)
        self.channel_type = 'Circular'
        self.D = D
        self.y = y if y != None else Symbol('y')
        if(type(self.y) is Symbol):
            self.calc_yn()
        if(self.Q is None):
            self.calc_flow()
        self.get_critical_parameters()
        

    def calc_properties(self):
        try:
            # Handle symbolic yn calculations first
            if self.yn is not None:
                try:
                    # Calculate and evaluate theta_n
                    theta_n_expr = acos((1 - (2 * (self.yn / self.D)))) * 2
                    # Store the symbolic expression for later use
                    self.theta_n = theta_n_expr
                    
                    # Get evaluated versions for calculations
                    theta_n_val = float(theta_n_expr.evalf())
                    sin_theta_n_val = float(sin(theta_n_expr).evalf())
                    
                    # Calculate other properties
                    self.ac = ((theta_n_val - sin_theta_n_val) * float(self.D)**2) / 8
                    self.pc = (float(self.D) * theta_n_val) / 2
                    self.twc = float(sin(theta_n_expr / 2).evalf()) * float(self.D)
                except (ValueError, TypeError):
                    self.ac = self.pc = self.twc = None
            
            # Handle numeric y calculations
            if self.y is not None and type(self.y) is not Symbol:
                try:
                    # Ensure we're working with float values
                    y = float(self.y)
                    d = float(self.D)
                    
                    # Sanity check to prevent math domain errors
                    depth_ratio = max(0, min(1, y / d))
                    angle_arg = 1 - 2 * depth_ratio
                    
                    # Calculate theta
                    theta_expr = acos(angle_arg) * 2
                    theta_val = float(theta_expr.evalf())
                    self.theta = theta_val
                    
                    # Calculate sin(theta) separately to avoid type errors
                    sin_theta_val = float(sin(theta_expr).evalf())
                    
                    # Calculate area
                    self.a = ((theta_val - sin_theta_val) * d**2) / 8
                    
                    # Calculate perimeter and other properties
                    self.p = (d * theta_val) / 2
                    
                    # Calculate the rest with proper error checking
                    if self.p > 0:
                        self.rh = self.a / self.p
                    else:
                        self.rh = None
                    
                    # Calculate top width
                    sin_half_theta_val = float(sin(theta_expr / 2).evalf())
                    self.tw = sin_half_theta_val * d
                    
                    # Calculate hydraulic depth
                    if self.tw > 0:
                        self.dh = self.a / self.tw
                    else:
                        self.dh = None
                    
                    # Calculate section factor
                    self.zc = d * (y**1.5)
                except (ValueError, TypeError, ZeroDivisionError) as e:
                    # Clear all properties if any calculation fails
                    self.a = self.p = self.rh = self.tw = self.dh = self.zc = None
            
            # Call parent method only if essential values are set
            if hasattr(self, 'a') and self.a is not None:
                super().calc_properties()
        except Exception as e:
            # Handle any unexpected errors
            self.a = self.p = self.rh = self.tw = self.dh = self.zc = None
        

    def __str__(self):
        return f"\nChannel: {self.channel_type}\nDimensions: \n\tDiameter: {self.D}\n{super().__str__()}"
    
    def caudal_manning(self):
        self.calc_properties()
        self.Q = (self.a * self.rh**(2/3) * self.So**0.5) / self.n
        self.calc_properties()
