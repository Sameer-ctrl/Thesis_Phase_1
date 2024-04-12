from astropy.io import fits
from astropy.table import Table
from numpy import *
from scipy.interpolate import interp2d,interp1d
import pickle
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings("ignore")

qso='3c263'
z_abs=0.140756

hdu=fits.open(f'{qso}/z={z_abs}/logZ=-1/component_II_PI_nH_Z_col_density_param.fits')
data=Table(hdu[1].data)

mask=data['H']>=10**13.45

ions=['Si+2','C+3','O+5']

Z=data['log_Z'][mask]
nH=data['log_nH'][mask]

for i in ions:

    plt.scatter(nH,Z,c=log10(data[i][mask]))
    plt.title(i)
    plt.colorbar()
    plt.show()


