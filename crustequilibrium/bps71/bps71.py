import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
import sys
from myers66 import *

# Pick a pressure P
# Pick Z,A pairs
# Calculate Gibbs free energy for each Z,A
#   Pre-drip:
#       Find n_N, n_e from P (given Z,A)
#       Calc mu(Z,A,n_N,n_e,P)
#   Post-drip:
#       ????????
# Find Z,A that give minimum mu at this P
# Increase P
# Repeat
# Now we have equilibrium Z,A at every P (or rho, etc.)

############

# GLOBAL stuff

mass_table = np.genfromtxt('crustequilibrium/bps71/mass_table.txt',dtype=None,skip_header=1)

num_nuclides = len(mass_table)

Zs = np.zeros(num_nuclides)
As = np.zeros(num_nuclides)
skynet_masses = np.zeros(num_nuclides)
names = []

for i in range(num_nuclides):
    Zs[i] = mass_table[i][0]
    As[i] = mass_table[i][1]
    skynet_masses[i] = mass_table[i][3]
    names.append(mass_table[i][2])

m_n = skynet_masses[0]  # MeV/c^2
m_p = skynet_masses[1]  # MeV/c^2
m_e = 0.511   # MeV/c^2

# For Myers66 masses
amu = 931.44   # MeV/c^2
m_n = amu + 8.07144
m_H = amu + 7.28899
m_p = m_H - m_e

MeV_to_grams = 1.7826619e-27 # times by this to go from MeV/c^2 -> grams
erg_to_MeV = 6.2415091e5  #  times by this to go from erg to MeV
esquared = (4.8032e-10)**2.0  # e^2 in cgs units (statC^2)
  # Note: e is in statcoulomb, which is dyne^0.5 * cm
  # Therefore e^2 is dyne * cm^2, and when dividing by 'a' in W_L, becomes
  # dyne * cm, which is erg
hbar = 1.0545716e-27  # erg s
hbar_MeV = hbar * erg_to_MeV
c = 2.99792458e10  #cm/s
threepisquared = (3.0 * np.pi**2.0)





def f1(x1,Z,P):  # solve for n_N (x) by rootfinding P - P_e - P_L = 0
    x = 10.0**x1
    
    n_e = Z*x
    #P_e_nonrel = 0.4 * threepisquared**(2.0/3.0) * \
    #             hbar**2.0 / (2.0 * m_e * MeV_to_grams) * \
    #             n_e**(5.0/3.0)
    P_e_rel = 0.25 * threepisquared**(1.0/3.0) * \
              (hbar*c) * n_e**(4.0/3.0)
    

    W_L = -1.81962*Z*Z * esquared * (x/2.0)**(1.0/3.0)
    P_L = (1.0/3.0) * W_L * x           

    #print(P, P_e, P_L)
    #return P - P_e_rel - P_L
    return P - P_e_BPS(n_e) - P_L


def Skynet_mass(i):
    #global skynet_masses
    
    return skynet_masses[i]# - Z*m_e
        # SkyNet mass table includes the (atomic) electrons right??
        # Actually no idea... Doesn't look like it for some tests (he4,fe56)
        

def equilibrium_predrip(P):

    mu = 1e99
    min_mu = 1e99
    min_ZA = [0,0]
    winner = {}
        
    for i in range(num_nuclides):
        Z = Zs[i]
        A = As[i]
        name = names[i]
        #print(Z,A)

        if Z >= 26 and Z <= 44: # can cheat to speed things up by
                                # imposing upper lim on tested Zs
            # 1e6 to 1e14 g/cc corresponds to nucleus number
            # densities around 1e6-1e14/(A*m_p) g/cc ~ 1e28-1e36 nuclei/cc
            
            #mass = Skynet_mass(i)
            mass = A*amu + Myers66_mass(Z,A) - Z*m_e
            if True:
                if name == 'fe56':
                    mass = mass + -0.717
                if name == 'ni62':
                    mass = mass + -0.446
                if name == 'ni64':
                    mass = mass + 0.145
                if name == 'ni66':
                    mass = mass + 0.694
                if name == 'kr86':
                    mass = mass + 0.886
            
            x = optimize.brentq(f1,20,40,args=(Z,P))
            n_N = 10.0**x
            n_e = Z * n_N

            W_N = mass  # in MeV
            W_L = (-1.81962*Z*Z* esquared * (n_N/2.0)**(1.0/3.0)) * erg_to_MeV

            # mu_e = (u_e + P_e)/n_e
            mu_e_rel = (1.0 * threepisquared**(1.0/3.0) * hbar*c * \
                        n_e**(1.0/3.0)) * erg_to_MeV
            mu_e_BPS = (u_e_BPS(n_e) + P_e_BPS(n_e))*erg_to_MeV/n_e
            


            mu = W_N + (4.0/3.0)*W_L + Z*mu_e_BPS
            mu = mu / A  # MeV
            if False: #name == 'ni62' or name == 'fe56':
                # debug stuff, use condition here to check specific nuclei
                print(name)
                print('W_N, (4.0/3.0)*W_L, Z*mu_e_BPS, mu')
                print(W_N,(4.0/3.0)*W_L,Z*mu_e_BPS,mu)
           

        if (mu < min_mu):
            min_mu = mu
            winner["Z"] = Z
            winner["A"] = A
            winner["mass"] = mass
            winner["name"] = name
            winner["n_N"] = n_N
            winner["n_e"] = n_e
            winner["mu"] = mu
            winner["mu_e"] = mu_e_BPS
            winner["rho"] = n_N * mass * MeV_to_grams
            
    #end for loop
        
    return winner

#end equilibrium_predrip

def u_e_BPS(n_e):
    k = (2.0 * 1.5 * np.pi * np.pi * n_e)**(1.0/3.0)
        # Factor of 2 required to match BPS results, likely subtle part of
        # definition of k and electrons having 2 spin states
    t = hbar * k / (m_e*MeV_to_grams * c)

    prefac = (m_e*MeV_to_grams)**4.0 * c**5.0 / (8.0 * np.pi*np.pi * hbar**3.0)

    return prefac * ( (2.0 * t**3.0 + t)*(t**2.0 + 1.0)**0.5 -
                      np.log(t + (t**2.0 + 1.0)**0.5)
                      )

def P_e_BPS(n_e):
    k = (2.0 * 1.5 * np.pi * np.pi * n_e)**(1.0/3.0)
    t = hbar * k / (m_e*MeV_to_grams * c)

    prefac = (m_e*MeV_to_grams)**4.0 * c**5.0 / (np.pi*np.pi * hbar**3.0)

    dudt = prefac * t**2.0 * (t**2.0 + 1.0)**0.5
    n_dudn = (1.0/3.0) * t * dudt

    return n_dudn - u_e_BPS(n_e)



if __name__ == '__main__':
    Ps = np.logspace(22,30,num=100) # dyne cm^-2, or erg cm^-3
    #Ps = np.logspace(23.7,23.7,num=2)

    with open("bps71_out.txt", "w") as myfile:
        myfile.write("{:>15}{:>12}{:>5}{:>5}{:>8}\n".format(
            "P [erg cm^-3]","rho [g/cc]","Z","A","mu_e"))

    #Start before neutron drip
    neutron_drip = False

    for i in range(len(Ps)):
        
        print('current P: '+str(Ps[i]))
        
        if neutron_drip == False:
            winner = equilibrium_predrip(Ps[i]) # This is where the work is done
            if winner["mu"] > m_n:
                neutron_drip = True
                print("Switching to neutron drip")
        if neutron_drip == True:
            print("no neutron regime yet")
            break
            
        if neutron_drip == False:
            print('** '+winner["name"]+'   ('+str(winner["Z"])+','+str(winner["A"])+
                  ')   rho = {:.2E} g/cc'.format(winner["rho"]))
            #print(winner["mu_e"])
            with open("bps71_out.txt", "a") as myfile:
                myfile.write("{:15.2E}{:12.2E}{:5.0f}{:5.0f}{:8.2f}\n".format(
                    Ps[i],winner["rho"],winner["Z"],winner["A"],winner["mu_e"]))

    # end loop over Ps




























####
####sys.exit()
##### All this is still WIP for now. Use with more robust version
####
####
####def get_mass(Z,A):
####    # This is for situations where you need to pull the mass
####    # for any Z,A pair, if not going through the entire
####    # list in order.
####    
####    Zindeces = [i for i,x in enumerate(Zs) if x == Z]
####    Aindeces = [i for i,x in enumerate(As) if x == A]
####
####    i = list(set(Zindeces).intersection(Aindeces))[0]
####
####    return(masses[i])
#####end get_mass
####
####
####
######def mass_MB77(Z,A):
######    N = A-Z
######    x = Z/A
######
######
######    k_0 = 1.34  # fm^-1
######    w_0 = 15.47  # MeV
######    s = 27.12  # MeV
######    K = 267.8  # MeV
######    sigma = 17.64  # MeV
######    delta = 0.2196  # fm
######
######    # Bulk energy
######    alpha = 1-2*x
######    W_onehalf = 
######
######
######    W_N = ((1-x)*m_n + x*m_p + W)*A + W_surf + W_Coul
####
