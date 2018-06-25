from crustequilibrium.bps71.bps71 import *


Ps = np.logspace(22,30,num=100) # dyne cm^-2, or erg cm^-3
#Ps = np.logspace(23.7,23.7,num=2)

with open("output/bps71_out.txt", "w") as myfile:
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
