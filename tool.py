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
    # PV_stats.print_stats()
    PV = Prism(PV_stats.diameter, PV_stats.diameter, PV_stats.length + PV_stats.diameter, 'PV', PV_stats.tot_weight_w_bat, PV_stats.tot_PV_disp_vol)
    return PV

def get_ratio_rel2width(params):
    w, l, h, hn, ht = params
    return [l/w, h/w, hn/w, ht/w]

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
    sys.exit(1)

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
        packer = Packer(f_params.tolist(), comps)
        res, runtime = packer.solve()
        vol = packer.compute_vol()
        print(f'Check result:   {res}')
        print(f'Fairing volume: {vol * 1e-9} m^3')
        print(f'Solve time:     {runtime} s')
        print()
        if res == unsat:
            f_params *= upscale
            record.append(packer)
            continue
        corners_outside_fairing = packer.corners_outside_fairing()
        record.append(packer)
        if corners_outside_fairing > 0:
            f_params *= upscale
        else: break

    print('*************************')
    print('       Downscaling       ')
    print('*************************')
    while True:
        comps[-1] = get_PV(f_params.tolist(), 892 * 1.2)
        packer = Packer(f_params.tolist(), comps)
        res, runtime = packer.solve()
        vol = packer.compute_vol()
        print(f'Check result:   {res}')
        print(f'Fairing volume: {vol * 1e-9} m^3')
        print(f'Solve time:     {runtime} s')
        print()
        if res == unsat:
            break
        corners_outside_fairing = packer.corners_outside_fairing()
        if packer.corners_outside_fairing() == 0:
            record.append(packer)
            f_params *= downscale
        else: break
    packer = record[-1]
    packer.plot()

    ratio_ub = 1.05
    ratio_lb = 0.95

    print('*************************')
    print('       Refinement        ')
    print('*************************')
    f_params = packer.fairing_params
    f_params_lb = np.array(f_params) * downscale
    vol = packer.compute_vol()
    rl, rh, rhn, rht = get_ratio_rel2width(f_params)
    s = Solver()
    # w = packer.fairing_params[0]
    w, l, h, hn, ht = Real('w'), Real('l'), Real('h'), Real('hn'), Real('ht')
    s.add(And(
        w * rl * ratio_lb <= l,
        w * rl * ratio_ub >= l))
    s.add(And(
        w * rh * ratio_lb <= h,
        w * rh * ratio_ub >= h))
    s.add(And(
        w * rhn * ratio_lb <= hn,
        w * rhn * ratio_ub >= hn))
    s.add(And(
        w * rht * ratio_lb <= ht,
        w * rht * ratio_ub >= ht))
    s.add(w  > f_params_lb[0])
    s.add(l  > f_params_lb[1])
    s.add(h  > f_params_lb[2])
    s.add(hn > f_params_lb[3])
    s.add(ht > f_params_lb[4])
    pack_success = True
    while True:
        if pack_success:
            s.add(w * l * (h + (hn + ht) / 3) <= vol * 0.95)
        
        if s.check() == sat:
            model = [z3RatNumRef2Float(s.model()[x]) for x in [w, l, h, hn, ht]]
            print(f'model: {model}')
            # comps[-1] = get_PV(model, 892 * 1.2)
            packer = Packer(model, comps)
            res, runtime = packer.solve()
            vol = packer.compute_vol()
            print(f'Check result:   {res}')
            print(f'Fairing volume: {vol * 1e-9} m^3')
            print(f'Solve time:     {runtime} s')
            print()
            if res == sat:
                pack_success = packer.corners_outside_fairing() == 0
                if pack_success:
                    record.append(packer)
                else:
                    s.add(w > model[0] * 1.02)
                    s.add(l > model[1] * 1.02)
                    s.add(h > model[2] * 1.02)
                    s.add(hn > model[3] * 1.02)
                    s.add(ht > model[4] * 1.02)
            else: break
        else: break

    print('End')
    record[-1].plot()
