import matplotlib.pyplot as plt
import numpy as np



rho_BPS = [8.1e6,2.7e8,1.2e9,8.2e9,2.2e10,4.8e10,1.6e11,
           1.8e11,1.9e11,2.7e11,3.7e11,4.3e11]
Z_BPS = [26,28,28,34,32,30,28,26,42, 40, 38, 36]
A_BPS = [56,62,64,84,82,80,78,76,124,122,120,118]
mu_e_BPS = [0.95,2.6,4.2,7.7,10.6,13.6,20.0,20.2,20.5,22.9,25.2,26.2]


a = np.genfromtxt('output/bps71_out.txt',skip_header=1)

P = a[:,0]
rho = a[:,1]
Z = a[:,2]
A = a[:,3]
mu_e = a[:,4]



## Z vs rho ##

for i in range(len(Z_BPS)):
    if i == 0:
        x = [1e0,rho_BPS[i]]
    else:
        deltarho = ((float(Z_BPS[i-1])/float(A_BPS[i-1]))
                    /(float(Z_BPS[i])/float(A_BPS[i])))
        x = [rho_BPS[i-1]*deltarho,rho_BPS[i]]
    y = [Z_BPS[i],Z_BPS[i]]
    plt.plot(x,y,color='black',linewidth=4.0,zorder=1)
#plt.scatter(rho_BPS[-1],Z_BPS[-1],s=6,color='red',marker='s')


plt.scatter(rho,Z,s=9,marker='.',color='dodgerblue',zorder=99)


plt.xlabel('rho [g/cc]')
plt.xlim([1e6,1e12])
plt.xscale('log')

plt.ylabel('Z')
plt.ylim([20,50])

#plt.show()

fig = plt.gcf()
fig.set_size_inches(11.0,7.0)
fig.savefig('output/bps71_out1.png',bbox_inches='tight')
plt.close(fig)
#############################



## A vs rho ##

for i in range(len(A_BPS)):
    if i == 0:
        x = [1e0,rho_BPS[i]]
    else:
        deltarho = ((float(Z_BPS[i-1])/float(A_BPS[i-1]))
                    /(float(Z_BPS[i])/float(A_BPS[i])))
        x = [rho_BPS[i-1]*deltarho,rho_BPS[i]]
    y = [A_BPS[i],A_BPS[i]]
    plt.plot(x,y,color='black',linewidth=4,zorder=1)
#plt.scatter(rho_BPS[-1],Z_BPS[-1],s=6,color='red',marker='s')


plt.scatter(rho,A,s=9,marker='.',color='dodgerblue',zorder=99)


plt.xlabel('rho [g/cc]')
plt.xlim([1e6,1e12])
plt.xscale('log')

plt.ylabel('A')
plt.ylim([50,130])

#plt.show()

fig = plt.gcf()
fig.set_size_inches(11.0,7.0)
fig.savefig('output/bps71_out2.png',bbox_inches='tight')
plt.close(fig)
#############################



# mu_e vs rho #

plt.scatter(rho_BPS,mu_e_BPS,s=12,color='black')

plt.plot(rho,mu_e,color='dodgerblue')

plt.xlabel('rho [g/cc]')
plt.xlim([1e6,1e12])
plt.xscale('log')

plt.ylabel(r'$\mu_e$ [MeV]')
plt.ylim([0,30])

#plt.show()

fig = plt.gcf()
fig.set_size_inches(11.0,7.0)
fig.savefig('output/bps71_out3.png',bbox_inches='tight')
plt.close(fig)
#############################




# Y_e vs rho #

for i in range(len(Z_BPS)):
    if i == 0:
        x = [1e0,rho_BPS[i]]
    else:
        deltarho = ((float(Z_BPS[i-1])/float(A_BPS[i-1]))
                    /(float(Z_BPS[i])/float(A_BPS[i])))
        x = [rho_BPS[i-1]*deltarho,rho_BPS[i]]
    Ye = float(Z_BPS[i])/float(A_BPS[i])
    y = [Ye,Ye]
    plt.plot(x,y,color='black',linewidth=3,zorder=1)

plt.plot(rho,np.divide(Z,A),color='dodgerblue')

plt.xlabel('rho [g/cc]')
plt.xlim([1e6,1e12])
plt.xscale('log')

plt.ylabel(r'$Y_e$')
plt.ylim([0.25,0.5])

#plt.show()

fig = plt.gcf()
fig.set_size_inches(11.0,7.0)
fig.savefig('output/bps71_out4.png',bbox_inches='tight')
plt.close(fig)
#############################



# mu_e vs rho #

for i in range(len(Z_BPS)):
    if i == 0:
        x = [0e0,mu_e_BPS[i]]
    else:
        x = [mu_e_BPS[i-1],mu_e_BPS[i]]
    Ye = float(Z_BPS[i])/float(A_BPS[i])
    y = [Ye,Ye]
    plt.plot(x,y,color='black',linewidth=3,zorder=1)

plt.plot(mu_e,np.divide(Z,A),color='dodgerblue')

plt.xlabel(r'$\mu_e$ [MeV]')
plt.xlim([0,30])

plt.ylabel(r'$Y_e$')
plt.ylim([0.25,0.5])
#plt.yscale('log')

#plt.show()

fig = plt.gcf()
fig.set_size_inches(11.0,7.0)
fig.savefig('output/bps71_out5.png',bbox_inches='tight')
plt.close(fig)
#############################




names_vs_Zs = (
        ['n',0],
        ['h',1],
        ['he',2],
        ['li',3],
        ['be',4],
        ['b',5],
        ['c',6],
        ['n',7],
        ['o',8],
        ['f',9],
        ['ne',10],
        ['na',11],
        ['mg',12],
        ['al',13],
        ['si',14],
        ['p',15],
        ['s',16],
        ['cl',17],
        ['ar',18],
        ['k',19],
        ['ca',20],
        ['sc',21],
        ['ti',22],
        ['v',23],
        ['cr',24],
        ['mn',25],
        ['fe',26],
        ['co',27],
        ['ni',28],
        ['cu',29],
        ['zn',30],
        ['ga',31],
        ['ge',32],
        ['as',33],
        ['se',34],
        ['br',35],
        ['kr',36],
        ['rb',37],
        ['sr',38],
        ['y',39],
        ['zr',40],
        ['nb',41],
        ['mo',42],
        ['tc',43],
        ['ru',44],
        ['rh',45],
        ['pd',46],
        ['ag',47],
        ['cd',48],
        ['in',49],
        ['sn',50],
        ['sb',51],
        ['te',52],
        ['i',53],
        ['xe',54],
        ['cs',55],
        ['ba',56],
        ['la',57],
        ['ce',58],
        ['pr',59],
        ['nd',60],
        ['pm',61],
        ['sm',62],
        ['eu',63],
        ['gd',64],
        ['tb',65],
        ['dy',66],
        ['ho',67],
        ['er',68],
        ['tm',69],
        ['yb',70],
        ['lu',71],
        ['hf',72],
        ['ta',73],
        ['w',74],
        ['re',75],
        ['os',76],
        ['ir',77],
        ['pt',78],
        ['au',79],
        ['hg',80],
        ['tl',81],
        ['pb',82],
        ['bi',83],
        ['po',84],
        ['at',85],
        ['rn',86],
        ['fr',87],
        ['ra',88],
        ['ac',89],
        ['th',90],
        ['pa',91],
        ['u',92],
        ['np',93],
        ['pu',94],
        ['am',95],
        ['cm',96],
        ['bk',97],
        ['cf',98],
        ['es',99],
        ['fm',100],
        ['md',101],
        ['no',102],
        ['lr',103],
        ['rf',104],
        ['db',105],
        ['sg',106],
        ['bh',107],
        ['hs',108],
        ['mt',109],
        ['ds',110],
        ['rg',111],
        ['cn',112],
        ['nh',113],
        ['fl',114],
        ['mc',115],
        ['lv',116],
        ['ts',117],
        ['og',118]
    )
