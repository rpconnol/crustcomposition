import numpy as np
import scipy.optimize as opt
import sys

from ..constants.constants_def import *
from ..mb77.mb77_lib import *

'''
adfadfa
'''


# REGARDING NEUTRON DRIP REGIME:
# Into the three equations in (5), go P, A_cell, and an (A,Z) pair.
# With those equations one should be able to self consistently calculate
# n_N, n_n, V_N, n_e, etc. To simplify things, let's pick the one thing
# we're going to boil everything down into: n_n. This works as long as I
# don't treat separately the W_L component of W_coul in the mb77 mass calc.
# W_L needs n_N explicitly, which turns it into a two simultaneous equation
# solve for n_n and n_N, if I want to include that contribution when solving
# for k. But the contribution is probably small relatively to the other terms
# in W_coul, so let's leave it separate for now, which means everything can be
# put in terms of n_n.

# Note: ALL units in MeV and fm for now. Densities are fm^-3, pressures are
# MeV/fm^3, k's are fm^-1, etc.


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

    mu_n = W_bulk(k_n,0) + k_n / 3.0 * dWdk(k_n,0)

    return n_n * mu_n - u_n(n_n)
#end P_n


def G_cell(P,A,Z,n_N,V_N,n_n):
    # Ensure pressure P is in MeV/fm^3
    
    n_e = Z * n_N
    
    return W_N_solve(A,Z,n_n) + W_L(Z,n_N,V_N) + \
           (u_e(n_e) + (1.0 - n_N*V_N)*u_n(n_n) + P)/n_N
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


def f_eq5_n_n(n_n_log10,P,A,Z,A_cell):
    n_n = 10.0**n_n_log10
    print(n_n)
    
    k_n = (1.5 * np.pi*np.pi * n_n)**onethird
    V_N = give_V_N(A,Z,k_n)

    n_N = ((A_cell - A)/n_n + V_N)**(-1.0)
    print(n_N)

    n_e = Z * n_N

    print(P,P_e(n_e),P_L(Z,n_N,V_N),P_n(n_n))
    return P - P_e(n_e) - P_L(Z,n_N,V_N) - P_n(n_n)
#def f_eq5


def f_eq5(n_N_log10,P,A,Z,A_cell):
    n_N = 10.0**n_N_log10
    print(n_N)

    n_e = Z * n_N
    
    k_n = 0.0
    for i in range(5):
        V_N = give_V_N(A,Z,k_n)
        
        n_n = (A_cell - A)/(n_N**(-1.0) - V_N)
        print(n_n)
        k_n = (1.5 * np.pi*np.pi * n_n)**onethird

    

    print(P,P_e(n_e),P_L(Z,n_N,V_N),P_n(n_n))
    return P - P_e(n_e) - P_L(Z,n_N,V_N) - P_n(n_n)
#def f_eq5


def innercrust_solve_n_n(P,A,Z,A_cell):
    # Calculates (n_N,V_N,n_n) for a given P, nucleus (Z,A), and A_cell

    try:
        n_n_log10 = opt.brentq(f_eq5,-10.0,-1.0,args=(P,A,Z,A_cell))
        n_n = 10.0**n_n_log10
        print(n_n)

        k_n = (1.5 * np.pi*np.pi * n_n)**onethird
        V_N = give_V_N(A,Z,k_n)

        n_N = ((A_cell - A)/n_n + V_N)**(-1.0)
        
        return [n_N,V_N,n_n]
    except ValueError:
        print('no n_n solution for ('+str(Z)+','+str(A)+')')
        return [0.0,0.0,0.0]
#end innercrust_solve


def innercrust_solve(P,A,Z,A_cell):
    # Calculates (n_N,V_N,n_n) for a given P, nucleus (Z,A), and A_cell

    n_N_log10 = opt.brentq(f_eq5,-20.0,-4.0,args=(P,A,Z,A_cell))
    n_N = 10.0**n_N_log10
    print(n_N)
        
    n_n = (A_cell - A)/(n_N**(-1.0) - V_N)
    k_n = (1.5 * np.pi*np.pi * n_n)**onethird
    
    V_N = give_V_N(A,Z,k_n)
    
    return [n_N,V_N,n_n]
#end innercrust_solve


# The big one
def equilibrium_ZA(P,A_cell,neutron_drip=False):
    Gibbs = 1e99
    min_Gibbs = 1e99
    min_ZA = [0,0]
        
    if neutron_drip == False:
        A = A_cell
        for Z in np.arange(10,np.ceil(A_cell/2)+1,dtype=int): # Test
            A = float(A)
            Z = float(Z)
            A_cell = float(A_cell)
            
            n_N = outercrust_solve(P,A,Z)
            V_N = 0.0
            n_n = 0.0

            Gibbs = G_cell(P,A,Z,n_N,V_N,n_n) / A_cell
            print(A,Z,Gibbs)
            if (Gibbs < min_Gibbs):
                min_Gibbs = Gibbs
                min_ZA = [Z,A]
                n_N_at_min = n_N
                n_n_at_min = n_n


    if neutron_drip == True:
        for A in np.arange(40,np.ceil(A_cell),dtype=int):
            for Z in np.arange(20,np.ceil(A/2),dtype=int):
                A = float(A)
                Z = float(Z)
                A_cell = float(A_cell)

                print(Z,A)
                [n_N,V_N,n_n] = innercrust_solve(P,A,Z,A_cell)

                
                if n_N > 0:
                    Gibbs = G_cell(P,A,Z,n_N,V_N,n_n) / A_cell
                    print(A,Z,Gibbs)
                    if (Gibbs < min_Gibbs):
                        min_Gibbs = Gibbs
                        min_ZA = [Z,A]
                        n_N_at_min = n_N
                        n_n_at_min = n_n
            

    # mass density in CGS units!
    rho_at_min = (n_N_at_min*1e39) * (W_N_solve(A,Z,n_n_at_min) * MeV_to_grams)
    
    return (min_ZA,rho_at_min)
#end equilibrium_ZA
    





        

