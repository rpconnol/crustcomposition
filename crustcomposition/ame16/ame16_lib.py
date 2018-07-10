import numpy as np


class AME16_TABLE:
    def __init__(self,user_file_path):
        self.ame16_table = self.load_mass16rc(user_file_path)
    #
    
    
    def load_mass16rc(self,file_path):
        ame16_data = np.genfromtxt(file_path,
                                   names=['A','Z','mass'],
                                   #dtype=['int','int','float'],
                                   skip_header=1)
        return ame16_data
    #


    def get_mass(self,Z,A):
        Zidx = np.where(self.ame16_table['Z'] == Z)
        Aidx = np.where(self.ame16_table['A'] == A)
        #print(Zidx)
        #print(Aidx)

        idx = [i for i in Zidx[0] if i in Aidx[0]]
        #print(idx)

        if idx == []:
            return -1
        else:
            return(float(self.ame16_table['mass'][idx]))
    #
    
## END CLASS ##




if __name__ == '__main__':
    mass_table = AME16_TABLE('../../data/mass16_rc.txt')

    print(mass_table.get_mass(2,4))
    print(mass_table.get_mass(1,14))
