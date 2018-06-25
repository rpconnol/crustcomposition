import numpy as np
import scipy.optimize as opt
import sys

from bbp71_constants import *
#from bbp71_solvers import *

'''
W_N(A,Z,k,k_n)
W_bulk(k,x)
w_surf(k,k_n,x)
W_prime(k,k_n,x)

dWdx(k,x)
d_Wprime_dx(k,k_n,x)
d_xwsurf_dx(k,k_n,x)
d_Wxequalszero_dnn(k)
d_Wthick_dnn(k,k_n,x)
d_wsurf_dnn(k,k_n,x)
dWdn(k,x)
d_Wprime_dn(k,k_n,x)
d_wsurf_dn(k,k_n,x)

little_v(k,k_n,x)
mu_e(k,k_n,x)
k_e(k,k_n,x)
mu_n_G(k,k_n,x)
P_G(k,k_n,x)
mu_n_N(k,k_n,x)
P_N(k,k_n,x)

A_solve(k,k_n,x)
mu_P_functions(k_in,x0)
kkn_solve(x,kkn_guess=[1.3,0.1])
'''


## W_N ##
def W_N(A,Z,k,k_n):
    A = float(A)
    Z = float(Z)

    x = Z/A

    #w_c,0 from w_c,0 Z^2 A^-1/3 term
    w_c0 = 0.6 * esquared_MeVfm * k / (9.0*np.pi/8.0)**onethird
    
    return ( (1.0-x)*m_n + x*m_p ) * A + \
           W_prime(k,k_n,x) * A + \
           w_surf(k,k_n,x) * A**twothirds + \
           w_c0 * Z**2.0 / A**onethird
## end W_N ##


## W(k,x) ##
def W_bulk(k,x):
    W_xequalszero = 19.74 * k**2.0 - k**3.0 * ((40.4 - 1.088*k**3.0)/
                                               (1.0 + 2.545*k))
    if x == 0:
        return W_xequalszero
    else:
        W_xequalsonehalf = -w_0 + K/(2.0*k_0**2.0) * (k-k_0)**2.0
        W_kin = 3.0 * (2.0**onethird * hbar_MeV * k)**2.0 / (10.0 * m_n_overc2) * \
                (x**(5.0/3.0) + (1.0-x)**(5.0/3.0))
        mu_p_0 = -1.0 * k**3.0 * ((218.0 + 277.0*k)/
                                  (1.0+8.57*k**2.0))
        mu_n_0 = W_xequalszero + k**2.0 / 3.0 * \
                 (2.13752*k**5.0 + 1.00787*k**4.0 + 7.73147*k**2.0 +
                  12.3132*k + 6.09539)/(k+0.392927)**2.0
        alpha = 1.0 - 2.0*x

        return (W_xequalsonehalf - 3.0*hbar_MeV**2.0*k**2.0/(10.0*m_n_overc2)) * \
               (1.0 - 3.0*alpha**4.0 + 2.0*alpha**6.0) \
               +\
               (s*(k/k_0)**2.0 - hbar_MeV**2.0*k**2.0/(6.0*m_n_overc2)) * \
               alpha**2.0 * (1.0-alpha**2.0)**2.0 \
               +\
               (W_xequalszero - 3.0*2.0**twothirds*hbar_MeV**2.0*k**2.0/(10.0*m_n_overc2)) * \
               (3.0*alpha**4.0 - 2.0*alpha**6.0) \
               +\
               0.25*(mu_p_0 - mu_n_0 + 2.0**twothirds*hbar_MeV**2.0*k**2.0/(2.0*m_n_overc2)) * \
               (alpha**4.0 - alpha**6.0) \
               +\
               W_kin
## end W_bulk ##


## w_surf ##
def w_surf(k,k_n,x):
    return sigma/w_0 * (W_bulk(k_n,0) - W_bulk(k,x)) * \
           (1.0 - (k_n/k)**3.0)**twothirds
## end w_surf ##


## W' = W(k,x) + W_thick + W_exch ##
def W_prime(k,k_n,x):
    
    # W_thick
    k_c = (k**3.0 - k_n**3.0)**onethird
    d = 0.74/k_c
    W_thick = -(4.0/9.0) * np.pi * x * esquared_MeVfm * d**2.0 * k**3.0 * x

    # W_exch
    W_exch = (-3.0/(4.0*np.pi)) * x * esquared_MeVfm * (2.0*x)**onethird * k

    return W_bulk(k,x) + W_thick + W_exch
## end W_prime ##





##### Derivatives #####


## d (W(k,x)) / dx ##              
def dWdx(k,x):
    alpha = 1.0 - 2.0*x
    
    W_xequalszero = 19.74 * k**2.0 - k**3.0 * ((40.4 - 1.088*k**3.0)/
                                               (1.0 + 2.545*k))
    W_xequalsonehalf = -w_0 + K/(2.0*k_0**2.0) * (k-k_0)**2.0
    mu_p_0 = -1.0 * k**3.0 * ((218.0 + 277.0*k)/
                              (1.0+8.57*k**2.0))
    mu_n_0 = W_xequalszero + k**2.0 / 3.0 * \
             (2.13752*k**5.0 + 1.00787*k**4.0 + 7.73147*k**2.0 +
              12.3132*k + 6.09539)/(k+0.392927)**2.0
    

    Wkx1 = W_xequalsonehalf - 3.0*hbar_MeV**2.0*k**2.0/(10.0*m_n_overc2)
    Wkx2 = s*(k/k_0)**2.0 - hbar_MeV**2.0*k**2.0/(6.0*m_n_overc2)
    Wkx3 = W_xequalszero - 3.0*2.0**twothirds*hbar_MeV**2.0*k**2.0/(10.0*m_n_overc2)
    Wkx4 = mu_p_0 - mu_n_0 + 2.0**twothirds*hbar_MeV**2.0*k**2.0/(2.0*m_n_overc2)
    dWkin_dx = (hbar_MeV * k)**2.0 / (2.0**onethird * m_n_overc2) * \
              (x**(2.0/3.0) - (1.0-x)**(2.0/3.0))

    return 24.0 * (Wkx1 - Wkx3) * (alpha**3.0 - alpha**5.0) \
           +\
           -4.0 * Wkx2 * (3.0*alpha**5.0 - 4.0*alpha**3.0 + alpha) \
           +\
           -1.0 * Wkx4 * (2.0*alpha**3.0 - 3.0*alpha**5.0) \
           +\
           dWkin_dx
## end Dwdx ##


## d (W') / dx ##
def d_Wprime_dx(k,k_n,x):
    k_c = (k**3.0 - k_n**3.0)**onethird
    d = 0.74/k_c
    return dWdx(k,x)  +\
           (-8.0/9.0) * np.pi * esquared_MeVfm * d**2.0 * k**3.0 * x  +\
           (-1.0/np.pi) * esquared_MeVfm * (2.0*x)**onethird * k
## end d_Wprime_dx ##


## d (x w_surf) / dx ##
def d_xwsurf_dx(k,k_n,x):
    return (sigma/w_0) * (1.0 - (k_n/k)**3.0)**twothirds * \
           (W_bulk(k_n,0) - W_bulk(k,x) - x*dWdx(k,x))
## end d_wsurf_dx ##


## d W(k_n,0) / dn_n ##
def d_Wxequalszero_dnn(k):
    return (0.5*np.pi*np.pi/k)*(2.13752*k**5.0 + 1.00787*k**4.0 +
                                7.73147*k**2.0 + 12.3132*k +
                                6.09539)/(k+0.392927)**2.0
## end d_Wxequalszero_dnn ##


## d W_thick / dn_n ##
def d_Wthick_dnn(k,k_n,x):
    return -0.73013 * np.pi**3.0 * esquared_MeVfm * k**3.0 * x**2.0 * \
           (k**3.0 - k_n**3.0)**(-3.0)
## end d_Wthick_dnn ##


## d w_surf / dn_n ##
def d_wsurf_dnn(k,k_n,x):
    return -1.0*sigma*np.pi*np.pi/(w_0 * k**3.0) * \
           (W_bulk(k_n,0) - W_bulk(k,x))/(1.0 - k_n**3.0 / k**3.0)*onethird
## end d_wsurf_dnn ##


## d W / dn ##
def dWdn(k,x):
    alpha = 1.0 - 2.0*x

    # derivatives with respect to k
    dW_xequalszero = (2.13752*k**6.0 + 1.00787*k**5.0 +
                      7.73147*k**3.0 + 12.3132*k**2.0 +
                      6.09539*k)/(k+0.392927)**2.0
    dW_xequalsonehalf = (K/k_0**2.0) * (k-k_0)
    dmu_p_0 = -1.0 * k**2.0 * (64.6441*k**3.0 + 25.4376*k**2.0 +
                               15.0861*k + 8.90463)/(k**2.0 + 0.116686)**2.0
    dmu_n_0 = (5.70006*k**9.0 + 9.63075*k**8.0 + 6.1163*k**7.0 +
               14.6148*k**6.0 + 33.8159*k**5.0 + 36.2336*k**4.0 +
               20.0152*k**3.0 + 5.57196*k**2.0 + 0.616292*k)*(k+0.392927)**(-5.0)
    dWkin_dk = 3.0 * (2.0**onethird * hbar_MeV)**2.0 * k / (5.0 * m_n_overc2) * \
              (x**(5.0/3.0) + (1.0-x)**(5.0/3.0))

    dWdk1 = (dW_xequalsonehalf - 0.6*hbar_MeV**2.0*k/m_n_overc2)
    dWdk2 = (s * 2.0*k/k_0**2.0 - onethird*hbar_MeV**2.0*k/m_n_overc2)
    dWdk3 = (dW_xequalszero - 0.6 * 2.0**twothirds*hbar_MeV**2.0*k/m_n_overc2)
    dWdk4 = dmu_p_0 - dmu_n_0 + 2.0**twothirds * hbar_MeV**2.0 * k / m_n_overc2

    dWdk = dWdk1 * (1.0 - 3.0*alpha**4.0 + 2.0*alpha**6.0) + \
           dWdk2 * alpha**2.0 * (1.0 - alpha**2.0)**2.0 + \
           dWdk3 * (3.0*alpha**4.0 - 2.0*alpha**6.0) + \
           0.25 * dWdk4 * (alpha**4.0 - alpha**6.0) + \
           dWkin_dk

    return 0.5 * np.pi*np.pi * k**(-2.0) * dWdk
## end dWdn ##


## d W' / dn ##
def d_Wprime_dn(k,k_n,x):
    return dWdn(k,x) - \
           0.47247 * np.pi * esquared_MeVfm * x**(4.0/3.0) * k**(-2.0) - \
           0.121689 * np.pi**3.0 * esquared_MeVfm * x**2.0 * \
           (k**3.0 - 3.0*k_n**3.0)/(k**3.0 - k_n**3.0)**(5.0/3.0)
## end d_Wprime_dn ##


## d w_surf / dn ##
def d_wsurf_dn(k,k_n,x):
    term1a = np.pi*np.pi * sigma * k_n**3.0 * (W_bulk(k_n,0) - W_bulk(k,x))
    term1b = w_0 * k**6.0 * (1.0 - k_n**3.0 / k**3.0)**onethird
    term2  = sigma/w_0 * (1.0 - k_n**3.0 / k**3.0)**twothirds * dWdn(k,x)

    return term1a/term1b - term2
## end d_wsurf_dn ##





##### Other #####


## v = (r_N/r_c)^3 = (0.5 / x) * (k_e / k)^3  [BBP 6.5] ##
def little_v(k,k_n,x):
    return (0.5/x) * (k_e(k,k_n,x) / k)**3.0
## end little_v ##


## mu_e from BBP 6.8 ##
def mu_e(k,k_n,x):
    A = A_solve(k,k_n,x)
    return (m_n - m_p) - d_Wprime_dx(k,k_n,x) - \
           x**(-1.0) * d_xwsurf_dx(k,k_n,x) * A**(-onethird)
## end mu_e ##

## k_e = mu_e / (hbar c) ##
def k_e(k,k_n,x):
    return mu_e(k,k_n,x)/(hbar_c_MeVfm)
## end k_e ##


## mu_n^(G) from BBP 6.13 ##
def mu_n_G(k,k_n,x):
    n = k**3.0 / (1.5 * np.pi**2.0)
    n_n = k_n**3.0 / (1.5 * np.pi**2.0)
    v = little_v(k,k_n,x)
    A = A_solve(k,k_n,x)

    return W_bulk(k_n,0) + n_n * d_Wxequalszero_dnn(k_n) + \
           n * v / (1.0-v) * (d_Wthick_dnn(k,k_n,x) +
                              A**(-onethird) * d_wsurf_dnn(k,k_n,x))
## end mu_N_G ##


## P^(G) from BBP 6.14 ##
def P_G(k,k_n,x):
    n_n = k_n**3.0 / (1.5 * np.pi**2.0)

    return n_n*(mu_n_G(k,k_n,x) - W_bulk(k_n,0))
## end P_G ##


## mu_n^(N) from BBP 6.16 ##
def mu_n_N(k,k_n,x):
    n = k**3.0 / (1.5 * np.pi**2.0)
    A = A_solve(k,k_n,x)

    return W_prime(k,k_n,x) + n * d_Wprime_dn(k,k_n,x) + \
           ((5.0/3.0)*w_surf(k,k_n,x) + n * d_wsurf_dn(k,k_n,x)) * A**(-onethird) + \
           x * (mu_e(k,k_n,x) - m_n + m_p)
## end mu_n_N ##


## P^(N) from BBP 6.18 ##
def P_N(k,k_n,x):
    n = k**3.0 / (1.5 * np.pi**2.0)
    v = little_v(k,k_n,x)
    A = A_solve(k,k_n,x)

    w_c0 = 0.6 * esquared_MeVfm * k / (9.0*np.pi/8.0)**onethird

    return n**2.0 * d_Wprime_dn(k,k_n,x) + n**2.0 * d_wsurf_dn(k,k_n,x) * A**(-onethird) + \
           onethird * n * w_c0 * x**2.0 * A**twothirds * (1.0 - v)
## end P_N ##









def A_solve(k,k_n,x):
    ### 1) Derivatives ###
    
    # In bbp71_terms now

    ### 2) lambda, nu, xi ###
    lambda_bbp = ((m_n - m_p) - d_Wprime_dx(k,k_n,x)) / \
                 (hbar_c_MeVfm * k * (2.0*x)**onethird)
    nu_bbp = d_xwsurf_dx(k,k_n,x) / \
             (hbar_c_MeVfm * k * 2.0**onethird * x**(4.0/3.0))
    xi_bbp = 2.5 * (np.pi/3.0)**onethird * w_surf(k,k_n,x) / \
             (x**2.0 * esquared_MeVfm * k)

    ### 3) a,b,c,d ###
    cubic_a = 2.0 - 3.0*lambda_bbp + lambda_bbp**3.0
    cubic_b = 3.0 * nu_bbp * (1.0 - lambda_bbp**2.0)
    cubic_c = 3.0 * nu_bbp**2.0 * lambda_bbp
    cubic_d = -1.0 * nu_bbp**3.0 - xi_bbp

    ### 4) Determinant, check < 0? ###
    cubic_determinant = 18.0*cubic_a*cubic_b*cubic_c*cubic_d - \
                        4.0 * cubic_b**3.0 * cubic_d + \
                        cubic_b**2.0 * cubic_c**2.0 - \
                        4.0 * cubic_a * cubic_c**3.0 - \
                        27.0 * cubic_a**2.0 * cubic_d**2.0
    #print("Determinant = "+str(cubic_determinant))

    ### 5) Delta0, Delta1, C ###
    cubic_delta0 = cubic_b**2.0 - 3.0 * cubic_a * cubic_c
    cubic_delta1 = 2.0 * cubic_b**3.0 - \
                   9.0*cubic_a*cubic_b*cubic_c + \
                   27.0 * cubic_a**2.0 * cubic_d
    cubic_bigC = -1.0*(-1.0*   # This is so we can take cube root of a negative (totally allowed)
                       (cubic_delta1 +
                        (-27.0 * cubic_a**2.0 * cubic_determinant)**0.5)/
                       2.0)**onethird
##    cubic_bigC_0 = (cubic_delta1 + (-27.0 * cubic_a**2.0 *
##                                    cubic_determinant)**0.5) / 2.0
##    if cubic_bigC_0 < 0:
##        cubic_bigC = -1.0*(-1.0*cubic_bigC_0)**onethird
##    else:
##        cubic_bigC = cubic_bigC_0**onethird

    cubic_bigC_minus = -1.0*(-1.0*
                             (cubic_delta1 -
                              (-27.0 * cubic_a**2.0 * cubic_determinant)**0.5)/
                             2.0)**onethird
##    cubic_bigC_minus_0 = (cubic_delta1 - (-27.0 * cubic_a**2.0 *
##                                          cubic_determinant)**0.5) / 2.0
##    if cubic_bigC_minus_0 < 0:
##        cubic_bigC_minus = -1.0*(-1.0*cubic_bigC_minus_0)**onethird
##    else:
##        cubic_bigC_minus = cubic_bigC_minus_0**onethird
    #print("C plus and C minus = ("+str(cubic_bigC)+"i,"+str(cubic_bigC_minus)+"i)")

    ### 6) A^1/3 solution from cubic ###
    A_onethird_root = (-1.0/(3.0*cubic_a)) * \
                      (cubic_b + cubic_bigC + cubic_delta0/cubic_bigC)

    ### 6b) A = A_onethird_root**3.0
    return A_onethird_root**3.0
## end A_solve ##


## mu_n_G - mu_n_N = 0 and P_G - P_N = 0 ##
def mu_P_functions(k_in,x0):
    k0 = k_in[0]
    k_n0 = k_in[1]
    if (k0 <= k_n0):  # if k = k_n, then k_c = 0 and things explode
        return [1e99,1e99]
    if (k0 <= 0) or (k_n0 <= 0):  # don't want negative k's
        return [1e99,1e99]
    if (k0 > 1.35):  # when k is too large, A_solve returns negative values
        return [1e99,1e99]

    f1 = mu_n_G(k0,k_n0,x0) - mu_n_N(k0,k_n0,x0)
    f2 = P_G(k0,k_n0,x0) - P_N(k0,k_n0,x0)

    return [f1,f2]
## end mu_P_functions ##


## Solve for k,k_n from x ##
def kkn_solve(x,kkn_guess=[1.3,0.1]):
    [k_out,k_n_out] = opt.fsolve(mu_P_functions,[kkn_guess[0],kkn_guess[1]],args=x)#,factor=0.1)

    return [k_out,k_n_out]
## end kkn_solve ##


## MANUAL k,k_n solve test (brute force) ##
def kkn_RC(x):
    k_arr = np.arange(0.01,2.00,0.01)
    k_n_arr = np.arange(0.01,2.00,0.01)

    k_sol = [0.0,0.0,1e99]
    for i in range(len(k_arr)):
        for j in range(len(k_n_arr)):
            k_out = k_arr[i]
            k_n_out = k_n_arr[j]
            [f1,f2] = mu_P_functions([k_out,k_n_out],x)
            if (abs(f1)+abs(f2)) < k_sol[2]:
                k_sol[0] = k_out
                k_sol[1] = k_n_out
                k_sol[2] = abs(f1)+abs(f2)

    return [k_sol[0], k_sol[1]]
## end kkn_RC ##










if __name__ == '__main__':
    print(hbar_c_MeV)
    
    # A_solve tests
    print('\n\nA(k,k_n,x)\n')

    print('Table 1 check')
    print('A(1.32,0.07,0.313) = '+str(A_solve(1.32,0.07,0.313)))
    print('A(1.32,0.12,0.310) = '+str(A_solve(1.32,0.12,0.310)))
    print('A(1.32,0.15,0.307) = '+str(A_solve(1.32,0.15,0.307)))
    print('A(1.30,0.36,0.280) = '+str(A_solve(1.30,0.36,0.280)))
    print('A(1.20,0.86,0.163) = '+str(A_solve(1.20,0.86,0.163)))
    print('A(1.30,1.25,0.057) = '+str(A_solve(1.30,1.25,0.057)))

    #print(kkn_RC(0.220))
    #sys.exit()

    # kkn_solve tests
    print('\nSolving for k,k_n from x')
    print('Table 1 check')
    print('x = 0.313:  '+str(kkn_solve(0.313)))
    print('x = 0.256:  '+str(kkn_solve(0.256)))
    print('x = 0.146:  '+str(kkn_solve(0.146)))
    print('x = 0.057:  '+str(kkn_solve(0.057)))

    # check to ensure P_G=P_N, mu_n_G=mu_n_N with the freshly solved k's
    print('\nP,mu crosschecks')
    x_check = 0.146
    print('x = '+str(x_check))
    [k_check,kn_check] = kkn_solve(x_check)
    print('solved (k,k_n) = ('+str(k_check)+', '+str(kn_check)+')')
    print('mu_n_G = '+str(mu_n_G(k_check,kn_check,x_check)))
    print('mu_n_N = '+str(mu_n_N(k_check,kn_check,x_check)))
    print('P_G = '+str(P_G(k_check,kn_check,x_check)))
    print('P_N = '+str(P_N(k_check,kn_check,x_check)))
    




    
    print('\nrho = 4.66e11; [A,Z,k,k_n] = [127,40,1.32,0.07]')
    print("(W_prime,w_surf) = ")
    print(W_prime(1.32,0.07,40.0/127.0),w_surf(1.32,0.07,40.0/127.0))
    print('Table 3: [-12.1, 15.4]')






    


    print('\n\n\nCalling all functions\n')
    
    [k,k_n,x] = [1.32,0.07,0.313]
    
    W_N(127.0,40.0,k,k_n)
    W_bulk(k,x)
    w_surf(k,k_n,x)
    W_prime(k,k_n,x)
    
    dWdx(k,x)
    d_Wprime_dx(k,k_n,x)
    d_xwsurf_dx(k,k_n,x)
    d_Wxequalszero_dnn(k)
    d_Wthick_dnn(k,k_n,x)
    d_wsurf_dnn(k,k_n,x)
    dWdn(k,x)
    d_Wprime_dn(k,k_n,x)
    d_wsurf_dn(k,k_n,x)
    
    little_v(k,k_n,x)
    mu_e(k,k_n,x)
    k_e(k,k_n,x)
    mu_n_G(k,k_n,x)
    P_G(k,k_n,x)
    mu_n_N(k,k_n,x)
    P_N(k,k_n,x)
    
    A_solve(k,k_n,x)
    mu_P_functions([k,k_n],x)
    kkn_solve(x)
    
    print('All clear')
        

