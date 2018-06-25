import matplotlib.pyplot as plt
import numpy as np


#Myers & Swiatecki 1966 fit parameters
a_1 = 15.677 #MeV
a_2 = 18.56 #MeV
c_3 = 0.717 #MeV
r_0 = 1.2049 #fm
c_4 = 1.21129 #MeV
kappa = 1.79
bigC = 5.8 #MeV
littleC = 0.26
a_over_r0 = 0.27


magic_numbers = [0,2,8,14,28,50,82,126,184]
def Myers66_shell_F(x):
    if x > magic_numbers[-1]:
        print("Myers66_shell_F("+str(x)+
              ") is outside largest magic number (184)")
        return 0.0
    for i in range(len(magic_numbers)-1):
        if x >= magic_numbers[i]:
            M_im1 = magic_numbers[i] # M_{i-1}, magic num below x
            M_i = magic_numbers[i+1] # M_i, next magic num above x
    #print("M_im1: "+str(M_im1))
    #print("M_i:   "+str(M_i))
    
    q = 0.6 * (M_i**(5.0/3.0) - M_im1**(5.0/3.0)) / (M_i - M_im1)
    #print("q:     "+str(q))

    # F(x), From equation below eq.5 in Myers66
    return q * (x - M_im1) - 0.6*(x**(5.0/3.0) - M_im1**(5.0/3.0))


def Myers66_mass(Z,A):
    Z = float(Z)
    A = float(A)
    m_n = 8.07144  # MeV, EXCESS
    m_H = 7.28899  # MeV, EXCESS
    #m_e = 0.511
    # defined at the top of the code:
    # a_1, a_2, c_3, c_4, r_0, kappa, bigC, littleC, a_over_r0
    N = A - Z
    c_1 = a_1 * ( 1.0 - kappa * ((N-Z)/A)**2.0 )
    c_2 = a_2 * ( 1.0 - kappa * ((N-Z)/A)**2.0 )
    #alpha_0_squared = 5.0 * a_over_r0**2.0 * (A**(-2.0/3.0))
    #alpha = a_2  # ???????

    # delta (in MeV)
    if int(A)%2 == 0:
        delta = -11.0 / A**0.5  # even-even
        if int(Z)%2 > 0:  # odd-odd
            delta = -1.0 * delta
    else:
        delta = 0.0 # even-odd or odd-even
    #print("delta = "+str(delta))

    # Shell term
    SNZ = bigC * ( (Myers66_shell_F(N) + Myers66_shell_F(Z)) /
                   (0.5 * A)**(2.0/3.0) -
                   littleC * A**(1.0/3.0)
                   )
    #print("S(N,Z) = "+str(SNZ))


    M = m_n*N + m_H*Z - \
        c_1 * A + \
        c_2 * A**(2.0/3.0) + \
        c_3 * Z**2.0 / A**(1.0/3.0) - \
        c_4 * Z**2.0 / A + \
        SNZ + \
        delta
    # note: omitted deformation stuff (alpha/gamma) in c_2, c_3, and SNZ terms
    #print("M = "+str(M)+" MeV")

    # note: M is the mass EXCESS. The rest mass would be
    #   A * 1 amu + M = A * 931.5 + M
    return M




# Begin tests

if __name__ == "__main__":

    if True:
        print("\n\nTEST Myers66_shell_F\n")

        print("\n-- Testing inbetweeners --\n")
        a = [1,5,12,22,30,55,99,150]
        for i in range(len(a)):
            print("Input: "+str(a[i]))
            result = Myers66_shell_F(a[i])
            print("Result: "+str(result))

        print("\n-- Testing flat q --\n")
        a = [51,55,59,67,80]
        for i in range(len(a)):
            print("Input: "+str(a[i]))
            result = Myers66_shell_F(a[i])
            print("Result: "+str(result))

        print("\n-- Testing magic numbers --\n")
        a = [2,50,126]
        for i in range(len(a)):
            print("Input: "+str(a[i]))
            result = Myers66_shell_F(a[i])
            print("Result: "+str(result))
                  
        print("\n-- Testing out of limits --\n")
        a = [200,300]
        for i in range(len(a)):
            print("Input: "+str(a[i]))
            result = Myers66_shell_F(a[i])
            print("Result: "+str(result))


    if True:
        print("\n\nTEST Myers66_mass\n")
        m_n = 939.51
        m_p = 938.22
        m_e = 0.511

        print("\nh1 (Z=2,A=4)")
        print(Myers66_mass(2,4))
        print("  (myers65 table: 7.483 MeV)")

        print("\ncr52 (Z=24,A=52)")
        print(Myers66_mass(24,52))
        print("  (myers65 table: -56.567 MeV)")

        print("\nfe56 (Z=26,A=56)")
        print(Myers66_mass(26,56))
        print("  (myers65 table: -59.888 MeV)")
        #print("  (wiki: 52103 MeV)")

        print("\npb208 (Z=82,A=208)")
        print(Myers66_mass(82,208))
        print("  (myers65 table: -20.722 MeV)")
        #print("  (wiki: 193729 MeV)")

        print("\n\n-- delta checks --")
        print("\nag107 (odd-even)")
        print(Myers66_mass(47,107))
        print("\nag108 (odd-odd)")
        print(Myers66_mass(47,108))
        
