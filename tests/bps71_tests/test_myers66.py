from crustequilibrium.bps71.myers66 import *


def run_tests():
    test_myers66_shell_F()
    test_myers66_mass()
###





def test_myers66_shell_F():
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


def test_myers66_mass():
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
    
