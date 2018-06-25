import numpy as np
import matplotlib.pyplot as plt

from bbp71_lib import *

x_tab1  = [0.313,0.310,0.307,0.304,0.302,0.299,0.295,0.291,0.286,0.280,
           0.273,0.266,0.256,0.246,0.234,0.220,0.202,0.186,0.163,0.146,
           0.121,0.096,0.081,0.067,0.057]
k_tab1  = [1.32,1.32,1.32,1.32,1.32,1.31,1.31,1.31,1.30,1.30,
           1.29,1.28,1.27,1.26,1.25,1.24,1.22,1.21,1.20,1.20,
           1.21,1.24,1.27,1.29,1.30]
kn_tab1 = [0.07,0.12,0.15,0.18,0.20,0.23,0.26,0.29,0.32,0.36,
           0.40,0.44,0.49,0.54,0.59,0.65,0.72,0.78,0.86,0.92,
           1.01,1.11,1.17,1.22,1.25]
A_tab1  = [127,130,134,137,140,144,149,154,161,170,
           181,193,211,232,257,296,354,421,548,683,
           990,1640,2500,4330,7840]


x = np.linspace(0.313,0.086,num=200)
k = np.zeros(len(x))
kn = np.zeros(len(x))
A = np.zeros(len(x))

for i in range(len(x)):
    if i <= 1:
        [k[i],kn[i]] = kkn_solve(x[i],kkn_guess=[1.3,0.1])
    else:
        #dkdx  = (k[i-1] - k[i-2]) / (x[i-1] - x[i-2])
        #dkndx = (kn[i-1] - kn[i-2])/(x[i-1] - x[i-2])
        #dx = x[i] - x[i-1]
        #[k[i],kn[i]] = kkn_solve(x[i],kkn_guess=[k[i-1]+dkdx*dx,kn[i-1]+dkndx*dx])
        [k[i],kn[i]] = kkn_solve(x[i],kkn_guess=[k[i-1],kn[i-1]])
    A[i] = A_solve(k[i],kn[i],x[i])



# k and k_n

plt.plot(x,k,color='k',label='$k$')
plt.plot(x,kn,color='blue',label='$k_n$')

plt.legend(loc=4,frameon=False,borderaxespad=2)

plt.scatter(x_tab1,k_tab1,color='k',s=20)
plt.scatter(x_tab1,kn_tab1,color='blue',s=20)
#plt.scatter(x_tab1,[y+0.005 for y in k_tab1],color='k',s=5)
#plt.scatter(x_tab1,[y+0.005 for y in kn_tab1],color='blue',s=5)
#plt.scatter(x_tab1,[y-0.005 for y in k_tab1],color='k',s=5)
#plt.scatter(x_tab1,[y-0.005 for y in kn_tab1],color='blue',s=5)

plt.plot([0.315,0.315],[0,1.4],color='gray',linestyle='--')
plt.annotate(s='neutron drip',xy=(0.313,0.7),xytext=(0.313,0.8),
             rotation='vertical',color='grey')
plt.plot([0.055,0.055],[0,1.4],color='gray',linestyle='--')
plt.annotate(s='core',xy=(0.055,0.7),xytext=(0.062,0.7),
             rotation='vertical',color='grey')

plt.xlim([0.32,0.05])
plt.ylim([0.0,1.4])

plt.xlabel('$x$')
plt.ylabel('$k$ [fm$^{-1}$]')

fig = plt.gcf()
fig.set_size_inches(8.0,5.0)
fig.savefig('bbp71_table1_a.png',bbox_inches='tight')
plt.close(fig)




# A

plt.plot(x,A,color='red',label='$A$')

#plt.legend(loc=4,frameon=False,borderaxespad=2)

plt.scatter(x_tab1,A_tab1,color='red',s=20)

plt.plot([0.315,0.315],[0,8000],color='gray',linestyle='--')
plt.annotate(s='neutron drip',xy=(0.313,0.7),xytext=(0.313,4500),
             rotation='vertical',color='grey')
plt.plot([0.055,0.055],[0,8000],color='gray',linestyle='--')
plt.annotate(s='core',xy=(0.055,0.7),xytext=(0.062,4000),
             rotation='vertical',color='grey')

plt.xlim([0.32,0.05])
plt.ylim([0,8000])

plt.xlabel('$x$')
plt.ylabel('$A$')

fig = plt.gcf()
fig.set_size_inches(8.0,5.0)
fig.savefig('bbp71_table1_b.png',bbox_inches='tight')
plt.close(fig)


