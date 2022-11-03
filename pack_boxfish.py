import sys
import time
import math
import numpy as np
from copy import deepcopy
from z3 import Solver, Real, Bool, And, Or, sat, unsat, Optimize

import stl_helper as SH
from helper import *
from plotter import Plotter
from prism import Prism

crush_depth = 1300
water_density = 1027
foam_density = 465


class PVCalculator:
    def __init__(self, length, diameter, energy_needed):
        self.diameter = diameter
        self.length = length
        diameter *= 1e-3
        length *= 1e-3
        g = 9.806
        youngs_modulus = 3.4e11
        yield_stress = 2.5e9
        poissons_ratio = 0.22
        material_density = 3.95
        crush_pressure = crush_depth * water_density * g
        
        self.energy_needed = energy_needed
        bat_weight_per_unit = 0.213
        vol_percentage = 0.67
        bat_vol = 0.0000987
        power_density_per_unit = 51.2
        power_density_per_mc = power_density_per_unit * vol_percentage / bat_vol / 1000
        # print(power_density_per_mc)
        num_bat_per_mc = vol_percentage / bat_vol
        num_bat_needed = energy_needed / power_density_per_mc * num_bat_per_mc
        # print(num_bat_needed)
        self.bat_mass = bat_weight_per_unit * num_bat_needed

        cylinder_thickness_elastic_failure = (crush_pressure / 2 / youngs_modulus * (1 - poissons_ratio ** 2)) ** (1/3) * diameter
        cylinder_thickness_yield_failure = 0.5 * (1 - math.sqrt(1 - (2 * crush_pressure / yield_stress))) * diameter
        cylinder_inner_diameter = diameter - 2 * max(cylinder_thickness_elastic_failure, cylinder_thickness_yield_failure)
        length = energy_needed / power_density_per_mc / (cylinder_inner_diameter / 2) ** 2 / math.pi
        self.length = length * 1e3
        cylinder_material_volume = math.pi * length * ((diameter / 2) ** 2 - (cylinder_inner_diameter / 2) ** 2)
        cylinder_weight = cylinder_material_volume * material_density * 1000

        k = crush_pressure / (0.69 * youngs_modulus) * math.sqrt(3 * (1 - poissons_ratio ** 2) / 4)
        endcap_thickness_elastic_failure = diameter * (k - 2 * math.sqrt(k)) / (k - 4)
        endcap_thickness_yield_failure = diameter / 2 - (diameter / 2) * (1 - 1.5 * (crush_pressure / yield_stress)) ** (1/3)
        endcap_inner_diameter = diameter - 2 * max(endcap_thickness_elastic_failure, endcap_thickness_yield_failure)
        endcap_material_volume = (4/3) * math.pi * ((diameter / 2) ** 3 - (endcap_inner_diameter / 2) ** 3)
        endcap_weight = endcap_material_volume * material_density * 1000

        self.tot_PV_weight = cylinder_weight + endcap_weight
        self.tot_PV_disp_vol = math.pi * ((diameter / 2) ** 2 * length + 4/3 * (diameter / 2) ** 3)
        # print(self.tot_PV_disp_vol)
        self.tot_buoy = self.tot_PV_disp_vol * water_density
        self.net_buoy = self.tot_buoy - self.tot_PV_weight

        self.tot_weight_w_bat = self.tot_PV_weight + self.bat_mass
        # print(self.bat_mass)

        # min_bat_diameter = 2 * math.sqrt(energy_needed / power_density_per_mc / length / math.pi)
        # print(min_bat_diameter)
        # assert cylinder_inner_diameter > min_bat_diameter, 'PV not large enough to pack battery'

    def print_stats(self):
        print('PV stats:')
        print(f'  diameter = {self.diameter} m')
        print(f'  length   = {self.length} m')
        print(f'  mass w/ bat = {self.tot_weight_w_bat} kg')
        print(f'  disp volume = {self.tot_PV_disp_vol} m^3')


class Packer:
    def __init__(self, fairing_params, comps, dyn_fairing=False, overlap=False):
        self.fairing_params_var = fairing_params 
        self.comps = deepcopy(comps)
        self.dyn_fairing = dyn_fairing
        self.overlap = overlap
        self.s = Solver()
        self.s = Optimize()
        self.res = None
        self.model = None
        self.fairing_params = None if dyn_fairing else fairing_params 
        self.fairing_vol = None
        self.fairing_disp_vol = None
        self.fairing_mass = None
        print(fairing_params)

        # NOTE: for static now only; dyanmic becomes SMTO problem
        if not self.dyn_fairing:
            f_w, f_l, f_h, f_hn, f_ht = self.fairing_params_var
            eps = 25
            thickness = 12.7
            fairing_density = 1450
            SH.gen_stl(f_h + eps, f_w + eps, f_l + eps, f_hn + eps, 0 + eps, f_ht + eps, 1)
            self.fairing_vol = SH.parse_volume()
            self.fairing_disp_vol = SH.parse_area() * thickness
            self.fairing_mass = self.fairing_disp_vol * fairing_density * 1e-9
            print(f'Fairing vol = {self.fairing_vol * 1e-9} m^3')
            print(f'Fairing area = {SH.parse_area() * 1e-6} m^2')
            print(f'Fairing disp = {self.fairing_disp_vol * 1e-9} m^3')
            print(f'Fairing mass = {self.fairing_mass} kg')

        self.encode()


    def encode(self):
        f_w, f_l, f_h, f_hn, f_ht = self.fairing_params_var

        # All corners in boxfish fairing
        for comp in self.comps:
            for i, coord in enumerate(comp.corners):
                x, y, z = coord
                in_fairing = Bool(f'{comp.name}_{i}')
                self.s.add(in_fairing == And(
                    f_ht * x - f_w / 2 * z <= 0,
                    f_ht * x + f_w / 2 * z >= 0,
                    f_ht * y - f_l / 2 * z <= 0,
                    f_ht * y + f_l / 2 * z >= 0,
                    x <= f_w / 2,
                    x >= -f_w / 2,
                    y <= f_l / 2,
                    y >= -f_l / 2,
                    f_hn * x + f_w / 2 * z <=  f_w / 2 * (f_h + f_ht + f_hn),
                    f_hn * x - f_w / 2 * z >= -f_w / 2 * (f_h + f_ht + f_hn),
                    f_hn * y + f_l / 2 * z <=  f_l / 2 * (f_h + f_ht + f_hn),
                    f_hn * y - f_l / 2 * z >= -f_l / 2 * (f_h + f_ht + f_hn)))
                self.s.add_soft(in_fairing)
                # self.s.add(in_fairing)


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
                no_overlap = Bool(f'o_{self.comps[i].name}_{self.comps[j].name}')
                self.s.add(no_overlap == Or(
                    corners_i[0][0] >= corners_j[-1][0],
                    corners_i[0][1] >= corners_j[-1][1],
                    corners_i[0][2] >= corners_j[-1][2],
                    corners_j[0][0] >= corners_i[-1][0],
                    corners_j[0][1] >= corners_i[-1][1],
                    corners_j[0][2] >= corners_i[-1][2],
                    ))
                # self.s.add_soft(overlap)
                self.s.add(no_overlap)

        # self.s.minimize(sum(overlap_vars))

        # Neutrally buoyant
        if not self.dyn_fairing:
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
            self.s.add(max_buoy == spare_volume * water_density + tot_comp_disp_vol * water_density)
            tot_weight = Real('tot_weight')
            self.s.add(tot_weight == tot_comp_mass + fairing_in_water_weight + spare_volume * foam_density)
            self.s.add(max_buoy >= tot_weight)

            residual = Real('residual')
            self.s.add(residual == max_buoy - tot_weight)
        
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
        PMT1 = next(c for c in self.comps if c.name == 'PMT1')
        PMT2 = next(c for c in self.comps if c.name == 'PMT2')
        self.s.add(PMT1.center[0] == -PMT2.center[0])
        self.s.add(PMT1.center[1] ==  PMT2.center[1])
        self.s.add(PMT1.center[2] ==  PMT2.center[2])

    def add_volume_constraint(self, vol):
        assert self.dyn_fairing
        f_w, f_l, f_h, f_hn, f_ht = self.fairing_params_var 
        self.s.add(f_w * f_l * ((f_hn + f_ht) / 3 + f_h) < vol)

    def solve(self):
        start = time.time()
        self.res = self.s.check()
        runtime = time.time() - start
        # print(self.res)
        if self.res == sat:
            self.model = self.s.model()
            model_dict = {}
            for x in self.model:
                model_dict[str(x)] = z3RatNumRef2Float(self.model[x])
            # print(model_dict)

            if self.dyn_fairing:
                self.fairing_params = [z3RatNumRef2Float(self.model[param_var]) for param_var in self.fairing_params_var]
            else:
                self.fairing_params = self.fairing_params_var
            print(self.fairing_params)

            self.corners_outside = 0
            for comp in self.comps:
                for i, corner in enumerate(comp.corners):
                    in_fairing = str(self.model[Bool(f'{comp.name}_{i}')]) == 'True'
                    if not in_fairing:
                        self.corners_outside += 1
            if self.corners_outside > 0:
                print(f'FAIL: {self.corners_outside} points outside fairing')
            else:
                print(f'Packed SUCCESSfully!')

        return self.res, runtime

    def compute_vol(self):
        assert self.fairing_params is not None and type(self.fairing_params[0]) is float
        f_w, f_l, f_h, f_hn, f_ht = self.fairing_params
        return f_w * f_l * ((f_hn + f_ht) / 3 + f_h)

    def corners_outside_fairing(self):
        assert self.res == sat
        return self.corners_outside


    def write_smt2(self):
        with open('tmp.smt2', 'w') as f:
            f.write(self.s.to_smt2())

    def plot(self):
        assert self.res == sat
        plotter = Plotter()

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

        plotter.pyramid(lower_rect, [0, 0, 0])
        plotter.pyramid(upper_rect, [0, 0, f_ht + f_h + f_hn])
        plotter.rect_prism(-f_w/2, -f_l/2, f_ht, f_w, f_l, f_h, color='b')
        
        # Plot components
        for comp in self.comps:
            comp.assign_loc(self.model)
            plotter.rect_prism(comp.x, comp.y, comp.z, comp.w, comp.l, comp.h)
            xs.extend([comp.x, comp.x + comp.w])
            ys.extend([comp.y, comp.y + comp.l])
            zs.extend([comp.z, comp.z + comp.h])

        for comp in self.comps:
            for i, corner in enumerate(comp.corners):
                in_fairing = str(self.model[Bool(f'{comp.name}_{i}')]) == 'True'
                if not in_fairing:
                    point = Prism.get_prism_corners(comp.x, comp.y, comp.z, comp.w, comp.l, comp.h)[i]
                    plotter.point(*point)

        low_mult, high_mult = 0.8, 1.2
        plotter.ax.set_xlim3d(low_mult * min(xs), high_mult * max(xs))
        plotter.ax.set_ylim3d(low_mult * min(ys), high_mult * max(ys))
        plotter.ax.set_zlim3d(low_mult * min(zs), high_mult * max(zs))
        plotter.ax.set_box_aspect((np.ptp(xs),np.ptp(ys),np.ptp(zs)))


        plotter.show()






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

    diameter = min(f_w, f_l) * 6/7
    PV_stats = PVCalculator(f_h - diameter, diameter, 892 * 1.2)
    PV_stats.print_stats()

    # PV = Prism(1400, 1400, 4700, 'PV', 5990+2175, 6.5167)
    PV = Prism(PV_stats.diameter, PV_stats.diameter, PV_stats.length + PV_stats.diameter, 'PV', PV_stats.tot_weight_w_bat, PV_stats.tot_PV_disp_vol)
    OBS = Prism(1000, 1000, 300, 'OBS', 91.3, 0.3)
    PMT1 = Prism(80, 80, 3000, 'PMT1', 75, 0.038)
    PMT2 = Prism(80, 80, 3000, 'PMT2', 75, 0.038)
    ALT = Prism(88, 88, 200, 'ALT', 60, 0)
    CTD = Prism(20, 160, 50, 'CTD', 60, 0)

    comps = [PV, OBS, PMT1, PMT2, ALT, CTD]
    # comps = [PV, OBS, PMT1, PMT2, ]

    packer = Packer([f_w, f_l, f_h, f_hn, f_ht], comps, dyn_fairing)
    res = sat
    iter = 1
    # while res == sat:
    print(f'Iteration {iter}:')
    res, runtime = packer.solve()
    if res == sat:
        vol = packer.compute_vol()
        print(f'Check result:   {res}')
        print(f'Fairing volume: {vol * 1e-9} m^3')
        print(f'Solve time:     {runtime} s')
        # packer.add_volume_constraint(vol)
        # packer.write_smt2()
        packer.plot()

    elif res == unsat:
        packer = Packer([f_w, f_l, f_h, f_hn, f_ht], comps, dyn_fairing, overlap=True)
        res, runtime = packer.solve()
        print(f'Check result:   {res}')
        # print(f'Fairing volume: {vol * 1e-9} m^3')
        print(f'Solve time:     {runtime} s')
        packer.plot()
