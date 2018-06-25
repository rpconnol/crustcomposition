from crustcomposition.bbp71.bbp71_lib import *

def run_tests():
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
