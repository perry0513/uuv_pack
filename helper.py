import z3
from z3 import Solver, Real, Bool, And, Or, sat, unsat, Optimize

def z3RatNumRef2Float(x):
    if type(x) == z3.RatNumRef:
        return x.numerator_as_long() / x.denominator_as_long()
    if type(x) == z3.BoolRef:
        return x
