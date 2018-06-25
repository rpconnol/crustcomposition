from crustcomposition.hz90.hz90_lib import *

def run_tests():
    neutron_drip = True

    test_pressure = 1e32
    test_pressure = test_pressure * erg_to_MeV / 1e39
    print(equilibrium_ZA(test_pressure,56,neutron_drip))

                                
##    if rho_at_min > 6e11 and neutron_drip == False:
##        neutron_drip = True
##        print("Switching to neutron drip mode")
