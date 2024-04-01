import matplotlib.pyplot as plt
from linetools.spectra.xspectrum1d import XSpectrum1D
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table,Column
from numpy import *
import os

qso='rxj0439'

file=f'Data/IGM_Danforth_Data/Spectra/{qso}_spec.fits'
file_systems=open(f'Data/IGM_Danforth_Data/Systems/{qso}_igm-systems.txt','r')
z=float(file_systems.readlines()[16].split(' ')[1])


hdu=fits.open(file)
data=Table(hdu[1].data)
wave=data['WAVE']
flux=data['FLUX']
err=data['ERR']


spec = XSpectrum1D.from_tuple((wave,flux,err))
spec.fit_continuum(kind='QSO', redshift=z)

cont=spec.co

tab=Table()

wave_col=Column(name='WAVE',data=wave)
flux_col=Column(name='FLUX',data=flux/cont)
err_col=Column(name='ERROR',data=err/cont)

tab.add_columns([wave_col,flux_col,err_col])
tab.write(f'Data/IGM_Danforth_Data/Cont_norm_spectra/{qso}_cont_norm.asc', format='ascii', overwrite=True)
tab.write(f'../VPfit/{qso}/{qso}_cont_norm.asc', format='ascii', overwrite=True)

os.system(f'cp _knots.jsn  ../VPfit/{qso}')

plt.step(wave,flux,label='spectrum')
plt.plot(wave,spec.co,label='continuum')
plt.legend()
plt.show()

# new_table.add_columns([wave_col,flux_col,err_col,cont_col,norm_flux_col,norm_err_col])
# new_table.write(f'{file[:-5]}_unbinned.fits', format='fits', overwrite=True)
# new_table.write(f'/home/sameer/Thesis_Phase_1/VPfit/PG0003+158_cont_norm.asc', format='ascii', overwrite=True)

# data['CONT_FLUX']=cont_col
# data['NORMALISED_FLUX']=norm_flux_col
# data['NORMALISED_ERROR']=norm_err_col
# data.add_columns([cont_col,norm_flux_col,norm_err_col])
# data.write(file,format='fits', overwrite='True')


# low_cont=cont*0.97
# cont_norm_flux=Column(name='FLUX',data=flux/low_cont)
# err_low=Column(name='ERROR',data=err/low_cont)

# tab=Table()
# tab.add_columns([wave,cont_norm_flux,err_low])
# tab.write(f'Data/{file[:-5]}_low_cont_norm.fits', format='fits', overwrite=True)
# tab.write(f'Data/{file[:-5]}_low_cont_norm.asc', format='ascii', overwrite=True)

# hi_cont=cont*1.03
# cont_norm_flux=Column(name='FLUX',data=flux/hi_cont)
# err_hi=Column(name='ERROR',data=err/hi_cont)

# tab=Table()
# tab.add_columns([wave,cont_norm_flux,err_hi])
# tab.write(f'Data/{file[:-5]}_high_cont_norm.fits', format='fits', overwrite=True)
# tab.write(f'Data/{file[:-5]}_high_cont_norm.asc', format='ascii', overwrite=True)
