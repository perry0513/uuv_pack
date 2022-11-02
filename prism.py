from z3 import Solver, Real, Bool, And, Or, sat, unsat, Optimize
from helper import *

'''
         _____ 
        /|   /|
       / |  / |
      /____/__|
      |  / |  /
      | /  | /
      |____|/
'''

class Prism:
    def __init__(self, w, l, h, name, mass, disp_vol):
        self.w = w
        self.l = l
        self.h = h
        self.name = name
        self.mass = mass
        self.disp_vol = disp_vol
        
        self.x_var = Real(f'{name}_x')
        self.y_var = Real(f'{name}_y')
        self.z_var = Real(f'{name}_z')
        self.x = None
        self.y = None
        self.z = None

        x = self.x_var
        y = self.y_var
        z = self.z_var

        self.corners = Prism.get_prism_corners(x, y, z, w, l, h)

        self.center = (x + w/2, y + l/2, z + h/2)
        
    @staticmethod
    def get_prism_corners(x, y, z, w, l, h): 
        xs = [x, x + w, x    , x    , x + w, x + w, x    , x + w]
        ys = [y, y    , y + l, y    , y + l, y    , y + l, y + l]
        zs = [z, z    , z    , z + h, z    , z + h, z + h, z + h]
        vs = list(zip(xs, ys, zs))
        return vs

    def assign_loc(self, model):
        self.x = z3RatNumRef2Float(model[self.x_var])
        self.y = z3RatNumRef2Float(model[self.y_var])
        self.z = z3RatNumRef2Float(model[self.z_var])
