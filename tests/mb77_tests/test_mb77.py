from crustequilibrium.mb77.mb77_lib import *

def run_tests():
    print(W_N_solve(40,20))
    print(W_N_solve(56,26))
    print('')
    print(give_me_BE(40,20))
    print(give_me_BE(56,26))
    print('')
    print(give_V_N(56,26,0.0))
    print(give_V_N(56,26,1.2))

    #for k in np.linspace(1.0,1.4,num=20):
    #    print(W_N(56.0,26.0,k,0.0))
