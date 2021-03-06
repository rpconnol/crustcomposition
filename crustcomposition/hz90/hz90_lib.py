import numpy as np
import scipy.optimize as opt
import sys

from ..constants.constants_def import *
from ..mb77.mb77_lib import *
from ..ame16.ame16_lib import *

'''

'''


mass_table = AME16_TABLE('data/mass16_rc.txt')




# Note: ALL units in MeV and fm for now. Densities are fm^-3, pressures are
# MeV/fm^3, k's are fm^-1, etc.


def get_mass(Z,A,n_n):
    mass = mass_table.get_mass(Z,A)
    if mass < 0:
        #print('No AME16 mass data for this nuclide, reverting to MB77')
        mass = W_N_solve(A,Z,n_n)
        #if mass < 0:
        #    print('-> No MB77 solution for this nuclide either!')

    return mass
#end get_mass


def get_VN(A,Z,k_n):
    if mass_table.get_mass(Z,A) > 0:
        #print('Nuclide exists in table, using fixed nuclear density')
        n_0 = k_0**3.0 / (1.5 * np.pi * np.pi) # saturation density from MB77's k_0
        return A/n_0
    else:
        #print('Nuclide not in table, calculating V_N from MB77')
        return give_V_N(A,Z,k_n)
#


def W_L(Z,n_N,V_N):
    # MB77 eq 5.10:
    # W_L  = -0.9 Z^2 e^2 / r_c * [ 1 - 0.333 * (r_N/r_c)^2 ]
    # with r_c = (4/3 pi n_N)^-1/3
    # W_L -> -0.9 * Z^2 * e^2 * (4/3 pi n_N)^1/3 * [ stuff ]
    #     -> -1.4508 * Z^2 * e^2 * n_N^1/3 * [ stuff ]
    # compare against BPS71:
    # W_L = -1.81962 Z^2 e^2 / a, n_N a^3 = 2
    #    -> -1.81962 Z^2 e^2 (n_N/2)^1/3
    #    -> -1.444 Z^2 e^2 n_N^1/3
    # so they're basically the same first order expression, the BBP71/MB77
    # version just has an additional correction for deeper crust when
    # r_N/r_c grows.

    r_N = (V_N / (4.0 / 3.0 * np.pi))**onethird
    r_c = (4.0 / 3.0 * np.pi * n_N)**(-1.0*onethird)

    return -1.45 * Z**2.0 * esquared_MeVfm * n_N**onethird * \
           (1.0 - onethird*(r_N/r_c)**2.0)
#end W_L


def P_L(Z,n_N,V_N=0.0):
    return (-1.0/3.0) * W_L(Z,n_N,V_N) * n_N
#end P_L


def P_e(n_e):
    return 0.25 * threepisquared**(1.0/3.0) * hbar_c_MeVfm * n_e**(4.0/3.0)
#end P_e


def u_e(n_e):
    # energy density of electrons (relativistic)
    return 0.75 * threepisquared**(1.0/3.0) * hbar_c_MeVfm * n_e**(4.0/3.0)
#end u_e


def u_n(n_n):
    # energy density of neutrons, including rest mass
    k_n = (1.5 * np.pi*np.pi * n_n)**onethird
    return W_bulk(k_n,0) * n_n + m_n * n_n
#end u_n


def P_n(n_n):
    # derived from eq.17 in BPS71,
    # P = n du/dn     - u
    #   = n d(n*W)/dn - (n*W)
    #   = n mu        - u
    k_n = (1.5 * np.pi*np.pi * n_n)**onethird

    #mu_n = W_bulk(k_n,0) + k_n / 3.0 * dWdk(k_n,0)
    # This expression for mu_n doesn't include the rest mass m_n,
    # because it's taken from the W_bulk calculation. However, u_n
    # DOES include the rest mass, so really du/dn for the pressure
    # SHOULD include rest mass if u_n does as well.

    #return n_n * mu_n - u_n(n_n)

    

    # So really I can just write the pressure as:
    # P = n du/dn - u
    #   = n (mu + m_n n_n) - u
    #   = n (mu + m_n n_n) - (W_bulk n_n + m_n n_n)
    #   = n * (1/3 k_n dWdk)

    return onethird * n_n * k_n * dWdk(k_n,0)
#end P_n


def G_cell(P,A,Z,n_N,V_N,n_n):
    # Ensure pressure P is in MeV/fm^3
    
    n_e = Z * n_N

    mass = get_mass(Z,A,n_n)
    if mass > 0:
        return mass + W_L(Z,n_N,V_N) + \
               (u_e(n_e) + (1.0 - n_N*V_N)*u_n(n_n) + P)/n_N
    else:
        return 1e99
#end G_cell


def f_eq3(n_N_log10,P,A,Z):
    n_N = 10.0**n_N_log10

    n_e = Z * n_N

    return P - P_e(n_e) - P_L(Z,n_N)
#def f_eq3


def outercrust_solve(P,A,Z):
    # Calculates n_N (and therefore n_e) for a given P and nucleus (Z,A)

    n_N_log10 = opt.brentq(f_eq3,-20,2,args=(P,A,Z))
    
    return 10.0**n_N_log10
#end outercrust_solve


##def f_eq5_n_n(n_n_log10,P,A,Z,A_cell):  # DEPRECATED
##    n_n = 10.0**n_n_log10
##    print(n_n)
##    
##    k_n = (1.5 * np.pi*np.pi * n_n)**onethird
##    V_N = give_V_N(A,Z,k_n)
##
##    n_N = ((A_cell - A)/n_n + V_N)**(-1.0)
##    print(n_N)
##
##    n_e = Z * n_N
##
##    print(P,P_e(n_e),P_L(Z,n_N,V_N),P_n(n_n))
##    return P - P_e(n_e) - P_L(Z,n_N,V_N) - P_n(n_n)
###def f_eq5_n_n


def f_eq5(n_N_log10,P,A,Z,A_cell):
    n_N = 10.0**n_N_log10
    #print(n_N)

    n_e = Z * n_N
    
    k_n = 0.0
    for i in range(3):
        V_N = get_VN(A,Z,k_n)
        
        n_n = (A_cell - A)/(n_N**(-1.0) - V_N)
        k_n = (1.5 * np.pi*np.pi * n_n)**onethird
#
    

    #print(P,P_e(n_e),P_L(Z,n_N,V_N),P_n(n_n))
    return P - P_e(n_e) - P_L(Z,n_N,V_N) - P_n(n_n)
#def f_eq5


##def innercrust_solve_n_n(P,A,Z,A_cell):   # DEPRECATED
##    # Calculates (n_N,V_N,n_n) for a given P, nucleus (Z,A), and A_cell
##
##    try:
##        n_n_log10 = opt.brentq(f_eq5,-10.0,-1.0,args=(P,A,Z,A_cell))
##        n_n = 10.0**n_n_log10
##        print(n_n)
##
##        k_n = (1.5 * np.pi*np.pi * n_n)**onethird
##        V_N = give_V_N(A,Z,k_n)
##
##        n_N = ((A_cell - A)/n_n + V_N)**(-1.0)
##        
##        return [n_N,V_N,n_n]
##    except ValueError:
##        print('no n_n solution for ('+str(Z)+','+str(A)+')')
##        return [0.0,0.0,0.0]
###end innercrust_solve


def innercrust_solve(P,A,Z,A_cell):
    # Calculates (n_N,V_N,n_n) for a given P, nucleus (Z,A), and A_cell

    try:
        n_N_log10 = opt.brentq(f_eq5,-20.0,-3.0,args=(P,A,Z,A_cell))
        n_N = 10.0**n_N_log10
        #print(n_N)
            
    ##    n_n = (A_cell - A)/(n_N**(-1.0) - V_N)
    ##    k_n = (1.5 * np.pi*np.pi * n_n)**onethird
    ##    
    ##    V_N = give_V_N(A,Z,k_n)

        k_n = 0.0
        for i in range(3):
            V_N = get_VN(A,Z,k_n)
            
            n_n = (A_cell - A)/(n_N**(-1.0) - V_N)
            #print(n_n)
            k_n = (1.5 * np.pi*np.pi * n_n)**onethird
        
        return [n_N,V_N,n_n]
    except ValueError:
        #print('no n_N solution for ('+str(Z)+','+str(A)+')')
        return [0.0,0.0,0.0]
#end innercrust_solve


# The big one
def equilibrium_ZA(P,A_cell,neutron_drip=False):
    Gibbs = 1e99
    min_Gibbs = 1e99
    [Z_sol,A_sol] = [0,0]
        
    if neutron_drip == False:
        A = A_cell
        for Z in np.arange(1,np.ceil(A_cell/2)+1,dtype=int): # Test
            A = float(A)
            Z = float(Z)
            A_cell = float(A_cell)

            n_N = outercrust_solve(P,A,Z)
            V_N = 0.0
            n_n = 0.0

            Gibbs = G_cell(P,A,Z,n_N,V_N,n_n) / A_cell
            #print(A,Z,Gibbs)
            if (Gibbs < min_Gibbs):
                min_Gibbs = Gibbs
                [Z_sol,A_sol] = [Z,A]
                n_N_sol = n_N
                n_n_sol = n_n


    if neutron_drip == True:
        for A in np.arange(1,np.ceil(A_cell)+1,dtype=int):
            for Z in np.arange(1,np.ceil(A/2)+1,dtype=int):
                A = float(A)
                Z = float(Z)
                A_cell = float(A_cell)

                #print(Z,A)
                [n_N,V_N,n_n] = innercrust_solve(P,A,Z,A_cell)
                #print(n_N,V_N,n_n)

                
                if n_N > 0:
                    Gibbs = G_cell(P,A,Z,n_N,V_N,n_n) / A_cell
                    #print(A,Z,Gibbs)
                    if (Gibbs < min_Gibbs):
                        min_Gibbs = Gibbs
                        [Z_sol,A_sol] = [Z,A]
                        n_N_sol = n_N
                        n_n_sol = n_n
            

    # mass density in CGS units!
    rho_sol = (n_N_sol*1e39) * \
                 (get_mass(Z_sol,A_sol,n_n_sol) * MeV_to_grams)

    # Y_n at solution
    Yn_sol = n_n_sol/(n_n_sol + n_N_sol*A_sol)

    # As a check, the number of free neutrons per nucleus should be
    # close to ~ A_cell-A (give or take, due to the nuclear volume)
    #print('free neutrons per nucleus: '+str(n_n_sol/n_N_sol))
    #print('Y_n = '+str(Yn_sol))

    return (Z_sol,A_sol,rho_sol,Yn_sol)
#end equilibrium_ZA
    


# the biggest one! Loop over the whole crust
def solve_crust(init_P,init_Acell,
                init_neutrondrip=False,P_stop=0.06,P_units='nuc'):

    # Start file
    with open("output/out.txt", "w") as myfile:
        myfile.write("Initial A_cell = {:0.0f}\n".format(init_Acell))
        myfile.write("{:>15}{:>12}{:>5}{:>5}{:>7}\n".format(
            "P [MeV fm^-3]","rho [g/cc]","Z","A","Y_n"))



    #initializing some things
    if P_units == 'nuc':
        P = init_P
    if P_units == 'cgs':
        P = init_P * erg_to_MeV / 1e39
        if P_stop > 1:  # nuclear P_stop should be something in 0.1-0.01 range
            P_stop = P_stop * erg_to_MeV / 1e39
    
    A_cell = float(init_Acell)

    neutron_drip = init_neutrondrip
    


    #Start looping over many pressures
    while P < P_stop:
        print('Attempting P = {:5.2E}'.format(P))
        (Z,A,rho,Y_n) = equilibrium_ZA(P,A_cell,neutron_drip)

        # Append to file
        with open("output/out.txt", "a") as myfile:
            myfile.write("{:15.2E}{:12.2E}{:5.0f}{:5.0f}{:7.2f}\n".format(
                P,rho,Z,A,Y_n))
        print('   Success, appended to file')


        if neutron_drip == False:
            # NEED A BETTER CONDITIONAL FOR TURNING ON DRIP MODE
            if rho > 6e11:
                neutron_drip = True
                print("Switching to neutron drip regime")

        
        P = 1.1 * P  # Maybe make this factor manually adjustable,
                     # or let user pick number of steps per decade or something
    #end while
    print('Stopping... At pressure limit')


    
#end solve_crust

