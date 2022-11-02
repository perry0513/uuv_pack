import sys
import time
import math
import numpy as np
from z3 import Solver, Real, Bool, And, Or, sat, unsat, Optimize

import stl_helper as SH
from helper import *
from plotter import Plotter
from prism import Prism
from pack_boxfish import Packer, PVCalculator

def get_PV(f_params, energy_needed):
    f_w, f_l, f_h, _, _= f_params
    diameter = min(f_w, f_l) * 6/7
    PV_stats = PVCalculator(f_h - diameter, diameter, energy_needed)
    PV_stats.print_stats()
    PV = Prism(PV_stats.diameter, PV_stats.diameter, PV_stats.length + PV_stats.diameter, 'PV', PV_stats.tot_weight_w_bat, PV_stats.tot_PV_disp_vol)
    return PV

def analyze(record):
    print('*************************')
    print('        Analysis         ')
    print('*************************')
    first = record[0]
    last = record[-1]
    if last.res == unsat:
        init_f_params = first.fairing_params
        print(f'> Initial fairing params: {init_f_params}')
        print(f'> Increase smallest dimension ({min(init_f_params)})')
    else:
        last.plot()
    pass

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
        f_params = np.array([f_w, f_l, f_h, f_hn, f_ht])



    # PV = Prism(1400, 1400, 4700, 'PV', 5990+2175, 6.5167)
    OBS = Prism(1000, 1000, 300, 'OBS', 91.3, 0.3)
    PMT1 = Prism(80, 80, 3000, 'PMT1', 75, 0.038)
    PMT2 = Prism(80, 80, 3000, 'PMT2', 75, 0.038)
    ALT = Prism(88, 88, 200, 'ALT', 60, 0)
    CTD = Prism(20, 160, 50, 'CTD', 60, 0)
    PV = get_PV(f_params.tolist(), 892 * 1.2)

    comps = [OBS, PMT1, PMT2, ALT, CTD, PV]
    # comps = [PV, OBS, PMT1, PMT2, ]

    upscale = 1.2
    downscale = 0.9
    param_ub = 10000

    res = unsat
    packer = None
    record = []
    print('*************************')
    print('        Upscaling        ')
    print('*************************')
    while True:
        if any(param > param_ub for param in f_params):
            print(f'Fairing parameter exceeds ub ({param_ub}), enter anaylsis...')
            analyze(record)
            break
        packer = Packer(f_params.tolist(), comps, dyn_fairing)
        res, runtime = packer.solve()
        if res == unsat:
            f_params *= upscale
            record.append(packer)
            continue
        corners_outside_fairing = packer.corners_outside_fairing()
        vol = packer.compute_vol()
        record.append(packer)
        print()
        print(f'Check result:   {res}')
        print(f'Fairing volume: {vol * 1e-9} m^3')
        print(f'Solve time:     {runtime} s')
        if corners_outside_fairing > 0:
            f_params *= upscale
        else: break

    print('*************************')
    print('       Downscaling       ')
    print('*************************')
    while True:
        packer = Packer(f_params.tolist(), comps, dyn_fairing)
        res, runtime = packer.solve()
        if res == unsat:
            break
        corners_outside_fairing = packer.corners_outside_fairing()
        vol = packer.compute_vol()
        print()
        print(f'Check result:   {res}')
        print(f'Fairing volume: {vol * 1e-9} m^3')
        print(f'Solve time:     {runtime} s')
        if packer.corners_outside_fairing() == 0:
            record.append(packer)
            f_params *= downscale
        else: break
    f_params = record[-1].fairing_params
    record[-1].plot()

        

