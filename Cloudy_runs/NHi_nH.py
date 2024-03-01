from astropy.table import Table
from numpy import *
from astropy.io import ascii,fits
import matplotlib.pyplot as plt

hdu=fits.open(f'NHi_nH/NHi_nH_col_density_param.fits')
data=Table(hdu[1].data)

nH=data['log_nH']
NHi=data['log_NHi']

ion='C+3'

N_Ovi=data['O+5']
N_ion=data[ion]

OVi_obs=13.63
SiIII_obs=12.39
Civ_obs=13.71

plt.scatter(log10(N_Ovi),log10(N_ion/N_Ovi),c=nH)
plt.colorbar()
plt.colorbar()

plt.scatter(OVi_obs,Civ_obs-OVi_obs,s=50,c='red')
plt.show()