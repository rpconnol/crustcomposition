import numpy as np



amu = 931.5   # MeV/c^2
m_n = 939.57  # MeV/c^2
m_p = 938.28  # MeV/c^2
m_e = 0.511   # MeV/c^2

MeV_to_grams = 1.7826619e-27 # times by this to go from MeV/c^2 -> grams
erg_to_MeV = 6.2415091e5  #  times by this to go from erg to MeV
fm_to_cm = 1e-13 # times by this to go from fm to cm
esquared = (4.8032e-10)**2.0  # e^2 in cgs units (statC^2)
  # Note: e is in statcoulomb, which is dyne^0.5 * cm
  # Therefore e^2 is dyne * cm^2, or erg * cm
esquared_MeVfm = esquared * erg_to_MeV / fm_to_cm
hbar = 1.0545716e-27  # erg s
hbar_MeV = hbar * erg_to_MeV  # MeV s
c = 2.99792458e10  # cm/s
c_fm = c/fm_to_cm  # fm/s
hbar_c = hbar * c  # erg cm
hbar_c_MeVfm = hbar_MeV * c_fm  # MeV fm
threepisquared = (3.0 * np.pi**2.0)
onethird = 1.0/3.0
twothirds = 2.0/3.0

# Fit parameters of MB77 (table 1)
k_0 = 1.34     # fm^-1
w_0 = 15.47    # MeV
s = 27.12      # MeV
K = 267.8      # MeV
sigma = 17.64  # MeV
delta = 0.2196 # fm
#

# More fit parameters from MB77 bulk energy section
# for W(k,0) and mu_p^(0)
c_sjoberg = [1.2974,
             15.0298,
             -15.2343,
             7.4663,
             23.1607,
             -142.37,
             38.6766,
             7.1178]
###





if __name__ == '__main__':
    print('esquared_MeVfm = '+str(esquared_MeVfm))
    print('hbar_MeV = '+str(hbar_MeV))
    print('c_fm = '+str(c_fm))
    print('hbar_c_MeVfm = '+str(hbar_c_MeVfm))
