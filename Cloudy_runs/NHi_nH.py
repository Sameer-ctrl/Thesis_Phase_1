from astropy.table import Table
from numpy import *
from astropy.io import ascii,fits
import matplotlib.pyplot as plt

hdu=fits.open(f'NHi_nH/NHi_nH_col_density_param.fits')
data=Table(hdu[1].data)

nH=data['log_nH']
NHi=data['log_NHi']

ion_Y='O+5'
ion_X='C+3'

obs_col_den={'C+3':13.71,'Si+2':12.39,'O+5':13.63}

def plot_grid(ion_X,ion_Y):

    N_ion_X=data[ion_X]
    N_ion_Y=data[ion_Y]



    # plt.scatter(log10(N_Ovi),log10(N_ion/N_Ovi),c=NHi)

    for n in unique(nH):
        mask=nH==n

        N_x_mask=N_ion_X[mask]
        N_y_mask=N_ion_Y[mask]
        plt.plot(log10(N_x_mask),log10(N_y_mask/N_x_mask),label=f'nH={n}')

    for n in unique(NHi):
        mask=NHi==n

        N_x_mask=N_ion_X[mask]
        N_y_mask=N_ion_Y[mask]
        plt.plot(log10(N_x_mask),log10(N_y_mask/N_x_mask),label=f'nH={n}')
        plt.scatter(obs_col_den[ion_X],obs_col_den[ion_Y]-obs_col_den[ion_Y],s=50,c='red')
        plt.xlabel(ion_X,fontsize=20)
        plt.ylabel(ion_Y,fontsize=20)
        plt.legend()

plt.figure(figsize=(40,20))

ax1=plt.subplot(121)
plot_grid('O+5','C+3')

plt.subplot(122,sharex=ax1, sharey=ax1)
plot_grid('O+5','Si+2')


plt.savefig('grid.png')
# plt.show()