amu = 931.4940954  # MeV
m_e = 0.5109989461 # MeV

Z_list = []
A_list = []
mass_list = []

with open('mass16.txt', 'r') as f:
    for i in range(39):
        line = f.readline()

    while True:
        # real data now
        line = f.readline()

        if line == '':
            #stop at the end of the file, empty line
            break

        Zstr = line[11:14]
        Astr = line[16:19]
        MEXstr = line[27:41]


        if ('#' not in MEXstr): # skip nuclei with non-experimental mass excess
            Z = float(Zstr)
            A = float(Astr)
            MEX = float(MEXstr) # keV
            #print(Zstr)
            #print(Z)
            #print(Astr)
            #print(A)
            #print(MEXstr)
            #print(MEX)
            

            atomic_mass = A*amu + MEX/1e3  # convert MEX from keV to MeV

            # total electron binding energy
            # approximate formula from Lunney, Pearson, and Thibault
            # (referenced via eq.2 in the AME2016-b paper)
            B_e = 14.4381 * Z**2.39 + 1.55468e-6 * Z**5.35  # in eV

            # Nuclear mass from atomic mass (AME2016-b eq.1)
            mass = atomic_mass - Z*m_e + B_e/1e6    # in MeV


            Z_list.append(Z)
            A_list.append(A)
            mass_list.append(mass)
            
# end looping over mass16.txt



with open("mass16_rc.txt", "w") as myfile:
    myfile.write("{:>3}{:>4}{:>15}\n".format(
        "A","Z","Nuc Mass (MeV)"))

    for i in range(len(A_list)):
        myfile.write("{:3.0f}{:4.0f}{:15.3f}\n".format(
            A_list[i],Z_list[i],mass_list[i]))

