import sys
import time
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
from z3 import Solver, Real, And, Or, sat, unsat

import helper as H

water_density = 1027
foam_density = 465

def embed():
    import IPython; IPython.embed()

def z3RatNumRef2Float(x):
    return x.numerator_as_long() / x.denominator_as_long()


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

        p1 = (x, y, z)
        p2 = (x + w, y, z)
        p3 = (x, y + l, z)
        p4 = (x, y, z + h)
        p5 = (x + w, y + l, z)
        p6 = (x + w, y, z + h)
        p7 = (x, y + l, z + h)
        p8 = (x + w, y + l, z + h)

        self.corners = [p1, p2, p3, p4, p5, p6, p7, p8]
        self.center = (self.x_var + w/2, self.y_var + l/2, self.z_var + h/2)

    def assign_loc(self, model):
        self.x = z3RatNumRef2Float(model[self.x_var])
        self.y = z3RatNumRef2Float(model[self.y_var])
        self.z = z3RatNumRef2Float(model[self.z_var])


class Packer:
    def __init__(self, fairing_params, comps, dyn_fairing):
        self.fairing_params_var = fairing_params 
        self.comps = comps
        self.dyn_fairing = dyn_fairing
        self.s = Solver()
        self.res = None
        self.model = None
        self.fairing_params = None
        self.plotter = Plotter()
        self.fairing_vol = None
        self.fairing_disp_vol = None
        self.fairing_mass = None

        # NOTE: for static now only; dyanmic becomes SMTO problem
        if not dyn_fairing:
            f_w, f_l, f_h, f_hn, f_ht = self.fairing_params_var
            eps = 25
            thickness = 12.7
            fairing_density = 1450
            H.gen_stl(f_h + eps, f_w + eps, f_l + eps, f_hn + eps, 0 + eps, f_ht + eps, 1)
            self.fairing_vol = H.parse_volume()
            self.fairing_disp_vol = H.parse_area() * thickness
            self.fairing_mass = self.fairing_disp_vol * fairing_density * 1e-9
            print(f'Fairing vol = {self.fairing_vol * 1e-9} m^3')
            print(f'Fairing area = {H.parse_area() * 1e-6} m^2')
            print(f'Fairing disp = {self.fairing_disp_vol * 1e-9} m^3')
            print(f'Fairing mass = {self.fairing_mass} kg')

        self.encode()


    def encode(self):
        f_w, f_l, f_h, f_hn, f_ht = self.fairing_params_var

        # All corners in boxfish fairing
        for comp in self.comps:
            for coord in comp.corners:
                x, y, z = coord
                self.s.add(f_ht * x - f_w / 2 * z <= 0)
                self.s.add(f_ht * x + f_w / 2 * z >= 0) 
                self.s.add(f_ht * y - f_l / 2 * z <= 0)
                self.s.add(f_ht * y + f_l / 2 * z >= 0)
                self.s.add(x <= f_w / 2)
                self.s.add(x >= -f_w / 2)
                self.s.add(y <= f_l / 2)
                self.s.add(y >= -f_l / 2)
                self.s.add(f_hn * x + f_w / 2 * z <=  f_w / 2 * (f_h + f_ht + f_hn))
                self.s.add(f_hn * x - f_w / 2 * z >= -f_w / 2 * (f_h + f_ht + f_hn))
                self.s.add(f_hn * y + f_l / 2 * z <=  f_l / 2 * (f_h + f_ht + f_hn))
                self.s.add(f_hn * y - f_l / 2 * z >= -f_l / 2 * (f_h + f_ht + f_hn))


        if self.dyn_fairing:
            self.s.add(f_w > 0)
            self.s.add(f_l > 0)
            self.s.add(f_h > 0)
            self.s.add(f_ht > 0)
            self.s.add(f_hn > 0)
            self.s.add(f_w < 10000)
            self.s.add(f_l < 10000)
            self.s.add(f_h < 10000)
            self.s.add(f_ht < 10000)
            self.s.add(f_hn < 10000)

        # No overlaps
        for i in range(len(self.comps)-1):
            for j in range(i+1, len(self.comps)):
                corners_i = self.comps[i].corners
                corners_j = self.comps[j].corners
                self.s.add(Or(
                    corners_i[0][0] >= corners_j[-1][0],
                    corners_i[0][1] >= corners_j[-1][1],
                    corners_i[0][2] >= corners_j[-1][2],
                    corners_j[0][0] >= corners_i[-1][0],
                    corners_j[0][1] >= corners_i[-1][1],
                    corners_j[0][2] >= corners_i[-1][2],
                    ))

        # Neutrally buoyant
        if not dyn_fairing:
            tot_comp_mass = Real('tot_comp_mass')
            self.s.add(tot_comp_mass == sum([c.mass for c in self.comps]))
            tot_comp_disp_vol = Real('tot_comp_disp_vol')
            self.s.add(tot_comp_disp_vol == sum([c.disp_vol for c in self.comps]))
            tot_comp_buoy = Real('tot_comp_disp_buoy')
            self.s.add(tot_comp_buoy == tot_comp_disp_vol * water_density)
            fairing_in_water_weight = Real('fairing_in_water_weight')
            self.s.add(fairing_in_water_weight == self.fairing_mass - self.fairing_disp_vol * water_density * 1e-9)

            spare_volume = Real('spare_volume')
            self.s.add(spare_volume == self.fairing_vol * 1e-9 - tot_comp_disp_vol)
            max_buoy = Real('max_buoy')
            self.s.add(max_buoy == spare_volume * (water_density - foam_density) + tot_comp_disp_vol * water_density)
            self.s.add(max_buoy >= tot_comp_mass + fairing_in_water_weight)

            residual = Real('residual')
            self.s.add(residual == max_buoy - tot_comp_mass - fairing_in_water_weight)
        
            # print(f'spare vol = {self.fairing_vol * 1e-9 - tot_comp_disp_vol}')
            # print(f'max_buoy = {(self.fairing_vol * 1e-9 - tot_comp_disp_vol) * (water_density - foam_density)}')
            # print(f'comp_mass + fairing_in_water = {tot_comp_mass + fairing_in_water_weight}')


            # Cb above cg
            # total_disp_volume = sum(
            tot_mass = tot_comp_mass + self.fairing_mass
            cm_x = sum([c.center[0] * c.mass / tot_mass for c in self.comps])
            cm_y = sum([c.center[1] * c.mass / tot_mass for c in self.comps])
            cm_z = sum([c.center[2] * c.mass / tot_mass for c in self.comps]) + self.fairing_mass / tot_mass * (f_ht + f_h + f_hn) / 2
            self.s.add(cm_x == 0)
            self.s.add(cm_y == 0)
            self.s.add(cm_z == f_ht + f_h / 2)

        # PMTs are symmetric
        self.s.add(self.comps[2].center[0] == -self.comps[3].center[0])
        self.s.add(self.comps[2].center[2] == self.comps[3].center[2])

    def add_volume_constraint(self, vol):
        assert self.dyn_fairing
        f_w, f_l, f_h, f_hn, f_ht = self.fairing_params_var 
        self.s.add(f_w * f_l * ((f_hn + f_ht) / 3 + f_h) < vol)

    def solve(self):
        start = time.time()
        self.res = self.s.check()
        runtime = time.time() - start
        print(self.res)
        if self.res == sat:
            self.model = self.s.model()
            model_dict = {}
            for x in self.model:
                model_dict[str(x)] = z3RatNumRef2Float(self.model[x])
            print(model_dict)

            if dyn_fairing:
                self.fairing_params = [z3RatNumRef2Float(self.model[param_var]) for param_var in self.fairing_params_var]
            else:
                self.fairing_params = self.fairing_params_var
            print(self.fairing_params)
        return self.res, runtime

    def compute_vol(self):
        assert self.fairing_params is not None and type(self.fairing_params[0]) is float
        f_w, f_l, f_h, f_hn, f_ht = self.fairing_params
        return f_w * f_l * ((f_hn + f_ht) / 3 + f_h)


    def write_smt2(self):
        with open('tmp.smt2', 'w') as f:
            f.write(self.s.to_smt2())

    def plot(self):
        assert self.res == sat

        xs, ys, zs = [], [], []

        f_w, f_l, f_h, f_hn, f_ht = self.fairing_params 
        
        # Plot fairing
        lower_rect = [
                    [-f_w/2, -f_l/2, f_ht],
                    [-f_w/2,  f_l/2, f_ht],
                    [ f_w/2,  f_l/2, f_ht],
                    [ f_w/2, -f_l/2, f_ht]]
        upper_rect = [
                    [-f_w/2, -f_l/2, f_ht + f_h],
                    [-f_w/2,  f_l/2, f_ht + f_h],
                    [ f_w/2,  f_l/2, f_ht + f_h],
                    [ f_w/2, -f_l/2, f_ht + f_h]]

        self.plotter.pyramid(lower_rect, [0, 0, 0])
        self.plotter.pyramid(upper_rect, [0, 0, f_ht + f_h + f_hn])
        self.plotter.rect_prism(-f_w/2, -f_l/2, f_ht, f_w, f_l, f_h, color='b')
        
        # Plot components
        for comp in self.comps:
            comp.assign_loc(self.model)
            self.plotter.rect_prism(comp.x, comp.y, comp.z, comp.w, comp.l, comp.h)
            xs.extend([comp.x, comp.x + comp.w])
            ys.extend([comp.y, comp.y + comp.l])
            zs.extend([comp.z, comp.z + comp.h])

        low_mult, high_mult = 0.8, 1.2
        self.plotter.ax.set_xlim3d(low_mult * min(xs), high_mult * max(xs))
        self.plotter.ax.set_ylim3d(low_mult * min(ys), high_mult * max(ys))
        self.plotter.ax.set_zlim3d(low_mult * min(zs), high_mult * max(zs))
        self.plotter.ax.set_box_aspect((np.ptp(xs),np.ptp(ys),np.ptp(zs)))

        self.plotter.show()

class Plotter:
    def __init__(self):
        fig = plt.figure()
        self.ax = fig.add_subplot(projection='3d')
        self.alpha = 0.2


    def rect_prism(self, x, y, z, w, l, h, color='r'):
        xs = [x, x + w, x    , x    , x + w, x + w, x    , x + w]
        ys = [y, y    , y + l, y    , y + l, y    , y + l, y + l]
        zs = [z, z    , z    , z + h, z    , z + h, z + h, z + h]
        vs = list(zip(xs, ys, zs))
        sides = [
                [vs[0], vs[1], vs[4], vs[2]],
                [vs[0], vs[1], vs[5], vs[3]],
                [vs[0], vs[2], vs[6], vs[3]],
                [vs[3], vs[5], vs[7], vs[6]],
                [vs[1], vs[4], vs[7], vs[5]],
                [vs[2], vs[4], vs[7], vs[6]]
                ]
        self.ax.add_collection3d(Poly3DCollection(sides, facecolors=color, linewidths=1, edgecolors=color, alpha=self.alpha))

    def pyramid(self, rect, peak, color='b'):
        sides = [
                [peak, rect[0], rect[1]],
                [peak, rect[1], rect[2]],
                [peak, rect[2], rect[3]],
                [peak, rect[3], rect[0]],
                [rect[0], rect[1], rect[2], rect[3]]
                ]
        self.ax.add_collection3d(Poly3DCollection(sides, facecolors=color, linewidths=1, edgecolors=color, alpha=self.alpha))

    def show(self):
        plt.show()





if __name__ == '__main__':
    dyn_fairing = len(sys.argv) != 6
    if dyn_fairing:
        f_w = Real('w')
        f_l = Real('l')
        f_h = Real('h')
        f_hn = Real('hn')
        f_ht = Real('ht')
    else:
        f_w, f_l, f_h, f_hn, f_ht = [float(arg) for arg in sys.argv[1:]]

    PV = Prism(1400, 1400, 4700, 'PV', 5990+2175, 6.5167)
    OBS = Prism(1000, 1000, 300, 'OBS', 91.3, 0.3)
    PMT1 = Prism(80, 80, 3000, 'PMT1', 75, 0.038)
    PMT2 = Prism(80, 80, 3000, 'PMT2', 75, 0.038)
    ALT = Prism(88, 88, 200, 'ALT', 60, 0)
    CTD = Prism(20, 160, 50, 'CTD', 60, 0)

    comps = [PV, OBS, PMT1, PMT2, ALT, CTD]
    # comps = [PV, OBS, ALT, CTD]

    packer = Packer([f_w, f_l, f_h, f_hn, f_ht], comps, dyn_fairing)
    res = sat
    iter = 1
    # while res == sat:
    print(f'Iteration {iter}:')
    res, runtime = packer.solve()
    vol = packer.compute_vol()
    print(f'Check result:   {res}')
    print(f'Fairing volume: {vol * 1e-9} m^3')
    print(f'Solve time:     {runtime} s')
    # packer.add_volume_constraint(vol)
    packer.plot()
