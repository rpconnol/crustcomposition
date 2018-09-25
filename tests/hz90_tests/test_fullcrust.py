from crustcomposition.hz90.hz90_lib import *

def run_tests():
    init_P = 1e26 # cgs
    init_Acell = 56 # fe56
    solve_crust(init_P,init_Acell,P_units='cgs')
