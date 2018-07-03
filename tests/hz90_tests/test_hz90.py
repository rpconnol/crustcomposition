from crustcomposition.hz90.hz90_lib import *

def run_tests():
    test_which_mass()
    test_pre_neutron_drip_singlepressure(1e29,56)
    test_neutron_drip_singlepressure(1e32,56)




def test_which_mass():
    print('Z,A = (26,56):')
    print(get_mass(26,56,0.0))

    print('\nZ,A = (26,80):')
    print(get_mass(26,80,0.0))



def test_pre_neutron_drip_singlepressure(test_pressure_cgs,A_cell):
    print('\nNo free neutrons')
    neutron_drip = False

    test_pressure = test_pressure_cgs * erg_to_MeV / 1e39
    print('Test pressure (MeV/fm): '+str(test_pressure))
    
    print('A_cell: '+str(A_cell))
        
    print(equilibrium_ZA(test_pressure,A_cell,neutron_drip))



def test_neutron_drip_singlepressure(test_pressure_cgs,A_cell):
    print('\nWith free neutrons')
    neutron_drip = True

    test_pressure = test_pressure_cgs * erg_to_MeV / 1e39
    print('Test pressure (MeV/fm): '+str(test_pressure))
    
    print('A_cell: '+str(A_cell))
        
    print(equilibrium_ZA(test_pressure,A_cell,neutron_drip))

                                
##    if rho_at_min > 6e11 and neutron_drip == False:
##        neutron_drip = True
##        print("Switching to neutron drip mode")
