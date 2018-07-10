import numpy as np
import scipy.optimize as opt
import sys

#from mb77_constants import *
from ..constants.constants_def import *

'''

'''

Wpair_in_WN = True



## W_N ##
def W_N(A,Z,k,k_n):
    A = float(A)
    Z = float(Z)

    x = Z/A

    if Wpair_in_WN == True:
        return ( (1.0-x)*m_n + x*m_p + W_bulk(k,x) ) * A + \
               W_surf(A,x,k,k_n) + W_coul(A,x,k,k_n) + W_pair(A,Z)
    else:
        return ( (1.0-x)*m_n + x*m_p + W_bulk(k,x) ) * A + \
               W_surf(A,x,k,k_n) + W_coul(A,x,k,k_n)
## end W_N ##


## W(k,x) ##
def W_bulk(k,x):
    W_xequalszero = c_sjoberg[0] * k      + \
                    c_sjoberg[1] * k**2.0 + \
                    c_sjoberg[2] * k**3.0 + \
                    c_sjoberg[3] * k**4.0
                    
    if x == 0:
        return W_xequalszero
    else:
        W_xequalsonehalf = 3.0*hbar_c_MeVfm**2.0*k**2.0/(10.0*m_n) * \
                           (1.0 - k/k_0)**3.0 - \
                           w_0 * (k/k_0)**3.0 * (1.0 + (1.0-k/k_0)*(9.0-6.0*k/k_0)) + \
                           0.5 * K * (1.0-k/k_0)**2.0 * (k/k_0)**3.0
        W_kin = 3.0 * (2.0**onethird * hbar_c_MeVfm * k)**2.0 / (10.0 * m_n) * \
                (x**(5.0/3.0) + (1.0-x)**(5.0/3.0))
        mu_p_0 = c_sjoberg[4] * k      + \
                 c_sjoberg[5] * k**2.0 + \
                 c_sjoberg[6] * k**3.0 + \
                 c_sjoberg[7] * k**4.0
        d_Wxequalszero_dk = c_sjoberg[0]                + \
                            2.0 * c_sjoberg[1] * k      + \
                            3.0 * c_sjoberg[2] * k**2.0 + \
                            4.0 * c_sjoberg[3] * k**3.0
        mu_n_0 = W_xequalszero + k / 3.0 * dWdk(k,0) #d_Wxequalszero_dk
        alpha = 1.0 - 2.0*x

        return (W_xequalsonehalf - 3.0*hbar_c_MeVfm**2.0*k**2.0/(10.0*m_n)) * \
               (1.0 - 3.0*alpha**4.0 + 2.0*alpha**6.0) \
               +\
               (s*(k/k_0)**2.0 - hbar_c_MeVfm**2.0*k**2.0/(6.0*m_n)) * \
               alpha**2.0 * (1.0-alpha**2.0)**2.0 \
               +\
               (W_xequalszero - 3.0*2.0**twothirds*hbar_c_MeVfm**2.0*k**2.0/(10.0*m_n)) * \
               (3.0*alpha**4.0 - 2.0*alpha**6.0) \
               +\
               0.25*(mu_p_0 - mu_n_0 + 2.0**twothirds*hbar_c_MeVfm**2.0*k**2.0/(2.0*m_n)) * \
               (alpha**4.0 - alpha**6.0) \
               +\
               W_kin
## end W_bulk ##


## W_surf ##
def W_surf(A,x,k,k_n=0.0):
    n_NM = k_0**3.0 / (1.5 * np.pi * np.pi)
    n_i  = k**3.0   / (1.5 * np.pi * np.pi)
    n_o  = k_n**3.0 / (1.5 * np.pi * np.pi)
    r_N  = (1.125 * np.pi)**onethird * A**onethird / k
    #print(A,A*x)
    #print(W_bulk(k_n,0),W_bulk(k,x))
    W_surf_0 = sigma * \
               (W_bulk(k_n,0) - W_bulk(k,x))**0.5 / w_0**0.5 * \
               (n_i - n_o)**2.0 / n_NM**2.0 * \
               k_0**2.0 / k**2.0 * \
               A**twothirds
    
    return W_surf_0 * (1.0 - 2.0 * delta / r_N)
## end W_surf ##


## W_coul ##
def W_coul(A,x,k,k_n=0.0):
    r_N = (1.125 * np.pi)**onethird * A**onethird / k
    d = 0.74/(k**3.0 - k_n**3.0)**onethird
    a1 = 5.0*np.pi*np.pi/6.0
    a2 = 18.031
    a3 = 1.3355

    W_coul_classic = 0.6 * (A*x)**2.0 ** esquared_MeVfm / r_N * \
                     (1.0 - a1 * (d/r_N)**2.0 + a2*(d/r_N)**3.0)
    W_coul_ex = -0.75 * (3.0/(2.0*np.pi))**twothirds * \
                (A*x)**(4.0/3.0) * esquared_MeVfm / r_N * \
                (1.0 - a3 * (d/r_N))

    #Intentionally leaving out the W_L contribution for now, as it's treated
    # separately in the HZ90 cell energy and pressure
    #n_N = (1e30 / A) * 1e-39
                # 1e30 baryons/cc, divided by A, convert to fm^-3
                # n_b = 1e30 cm^-3 is roughly ~2e6 g/cc of iron-56, from bps71
    #r_c = (1.0 / (1.333 * np.pi * n_N))**onethird
    #W_L = -0.9 * (A*x)**2.0 * esquared_MeVfm / r_c * \
    #      (1.0 - onethird*(r_N/r_c)**2.0)

    return W_coul_classic + W_coul_ex #+ W_L
## end W_coul ##


## d W / dk ##
def dWdk(k,x):
    d_Wxequalszero_dk = c_sjoberg[0]                + \
                        2.0 * c_sjoberg[1] * k      + \
                        3.0 * c_sjoberg[2] * k**2.0 + \
                        4.0 * c_sjoberg[3] * k**3.0
    
    if x == 0:
        return d_Wxequalszero_dk
    else:
        # derivatives with respect to k
        dW_xequalsonehalf = k_0**(-5.0) * (3.0*hbar_c_MeVfm**2.0*k**2.0/(10.0*m_n) *
                                           k_0**2.0 * (-3.0*k_0**2.0 + 6.0*k_0*k - 3.0*k**2.0) +
                                           k**2.0 * (k_0**2.0 * (1.5*K - 30.0*w_0) +
                                                     k_0 * k * (60.0*w_0 - 4.0*K) +
                                                     k**2.0 * (2.5*K - 30.0*w_0)))
        dmu_p_0 = c_sjoberg[4]      + \
                  2.0 * c_sjoberg[5] * k + \
                  3.0 * c_sjoberg[6] * k**2.0 + \
                  4.0 * c_sjoberg[7] * k**3.0
        dmu_n_0 = d_Wxequalszero_dk + onethird * (c_sjoberg[0] +
                                                  4.0 * c_sjoberg[1] * k +
                                                  9.0 * c_sjoberg[2] * k**2.0 +
                                                  16.0 * c_sjoberg[3] * k**3.0)
        dWkin_dk = 3.0 * (2.0**onethird * hbar_c_MeVfm)**2.0 * k / (5.0 * m_n) * \
                  (x**(5.0/3.0) + (1.0-x)**(5.0/3.0))

        dWdk1 = (dW_xequalsonehalf - 0.6*hbar_c_MeVfm**2.0*k/m_n)
        dWdk2 = (s * 2.0*k/k_0**2.0 - onethird*hbar_c_MeVfm**2.0*k/m_n)
        dWdk3 = (d_Wxequalszero_dk - 0.6 * 2.0**twothirds*hbar_c_MeVfm**2.0*k/m_n)
        dWdk4 = dmu_p_0 - dmu_n_0 + 2.0**twothirds * hbar_c_MeVfm**2.0 * k / m_n
        alpha = 1.0 - 2.0*x

        return dWdk1 * (1.0 - 3.0*alpha**4.0 + 2.0*alpha**6.0) + \
               dWdk2 * alpha**2.0 * (1.0 - alpha**2.0)**2.0 + \
               dWdk3 * (3.0*alpha**4.0 - 2.0*alpha**6.0) + \
               0.25 * dWdk4 * (alpha**4.0 - alpha**6.0) + \
               dWkin_dk
## end dWdk ##


def W_pair(A,Z):
    return -0.5 * ( (-1.0)**int(A-Z) + (-1.0)**int(Z) ) * 11.0/A**onethird
## end W_pair ##


def mb77_2p9(k,A,x,k_n):
    return k * dWdk(k,x) - (2.0 * W_surf(A,x,k,k_n) - W_coul(A,x,k,k_n))/A



def k_solve(A,Z,k_n):
    x = Z/A
    
    k = opt.brentq(mb77_2p9,max(1.2,k_n+0.01),1.4,args=(A,x,k_n))
    #print(k)
    
    return k




def W_N_solve(A,Z,n_n=0.0):
    Z = float(Z)
    A = float(A)

    k_n = (1.5 * np.pi*np.pi * n_n)**onethird

    try:
        k = k_solve(A,Z,k_n)
        
        if False:
            print("solved k  : "+str(k))
            print("n from k  : "+str(k**3.0 / (1.5 * np.pi * np.pi)))
            print("sum of p+n: "+str(Z*m_p+(A-Z)*m_n))
            print("W(k,x) * A: "+str(W_bulk(k,Z/A)*A))
            print("W_surf    : "+str(W_surf(A,Z/A,k,k_n)))
            print("W_coul    : "+str(W_coul(A,Z/A,k,k_n)))
            #print("W_pair    : "+str(-0.5 * ( (-1.0)**int(A-Z) + (-1.0)**int(Z) ) * 11.0/A**onethird))

        return W_N(A,Z,k,k_n)
    except ValueError:
        return -1



def give_me_BE(A,Z,k_n=0.0):
    Z = float(Z)
    A = float(A)
    
    mass = W_N_solve(A,Z)

    #W_pair = -0.5 * ( (-1.0)**int(A-Z) + (-1.0)**int(Z) ) * 11.0/A**onethird
    #mass = mass + W_pair

    BE = Z*m_p + (A-Z)*m_n - mass
    BE_per_nuc = BE/A

    return BE_per_nuc


def give_V_N(A,Z,k_n):
    Z = float(Z)
    A = float(A)

    k = k_solve(A,Z,k_n)

    # n = k^3 / (1.5 pi^2)
    # n = A/V_N -> V_N = A/n = 1.5 * pi^2 * A / k^3

    # returns V_N in fm^3
    return 1.5 * np.pi * np.pi * A / k**3.0
#end give_V_N

    
##def give_V_N(A,Z,k_n):
##    Z = float(Z)
##    A = float(A)
##
##    try:
##        k = k_solve(A,Z,k_n)
##
##        # n = k^3 / (1.5 pi^2)
##        # n = A/V_N -> V_N = A/n = 1.5 * pi^2 * A / k^3
##
##        # returns V_N in fm^3
##        return 1.5 * np.pi * np.pi * A / k**3.0
##    except ValueError:
##        print('ValueError in give_V_N')
##        return 0.0
###end give_V_N





