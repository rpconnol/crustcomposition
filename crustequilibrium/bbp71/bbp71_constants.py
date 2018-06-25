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

# m_n in MeV/c^2 with the c^2 divided through, specifically for use in
# [hbar^2 k^2 / m]  terms to cancel out the s^2 fm^-2 portion of the units
m_n_overc2 = m_n/(c_fm**2.0)  

# Fit parameters
w_0 = 16.5   # MeV
k_0 = 1.43   # fm^-1
K = 143.0    # MeV
s = 33.0     # MeV
sigma = 21.0 # MeV
# Note: because of k_0, all incoming k should be fm^-1
