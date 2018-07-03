import numpy as np


def load_mass16rc(file_path):
    ame16_data = np.genfromtxt(file_path,
                               names=['A','Z','mass'],
                               #dtype=['int','int','float'],
                               skip_header=1)
    return ame16_data
#


def mass_from_ZA(Z,A,ame16_data):
    Zidx = np.where(ame16_data['Z'] == Z)
    Aidx = np.where(ame16_data['A'] == A)
    #print(Zidx)
    #print(Aidx)

    idx = [i for i in Zidx[0] if i in Aidx[0]]
    #print(idx)

    if idx == []:
        return -1
    else:
        return(float(ame16_data['mass'][idx]))
#




if __name__ == '__main__':
    ame16_data = load_mass16rc('mass16_rc.txt')

    print(mass_from_ZA(2,4,ame16_data))
    print(mass_from_ZA(1,14,ame16_data))
