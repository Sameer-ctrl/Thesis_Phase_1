import matplotlib.pyplot as plt
from linetools.spectra.xspectrum1d import XSpectrum1D
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table,Column


file='Data/PG0003+158.fits'

hdu=fits.open(file)
data=Table(hdu[1].data)
wave=data['WAVE'][data['WAVE']>=1132.7]
flux=data['FLUX'][data['WAVE']>=1132.7]
err=data['ERROR'][data['WAVE']>=1132.7]

spec = XSpectrum1D.from_tuple((wave,flux,err))
spec.fit_continuum(kind='QSO', redshift=0.4509)


# plt.step(wave,flux,label='spectrum')
# plt.plot(wave,spec.co,label='continuum')
# plt.legend()
# plt.show()

cont_col=Column(name='CONT_FLUX',data=spec.co)

new_table=Table()

wave_col=Column(name='WAVE',data=wave)
flux_col=Column(name='FLUX',data=flux)
err_col=Column(name='ERROR',data=err)
norm_flux_col=Column(name='NORMALISED_FLUX',data=flux/spec.co)
norm_err_col=Column(name='NORMALIZED_ERROR',data=err/spec.co)

new_table.add_columns([wave_col,flux_col,err_col,cont_col,norm_flux_col,norm_err_col])
new_table.write(f'{file[:-5]}_unbinned.fits', format='fits', overwrite=True)
new_table.write(f'{file[:-5]}_unbinned.asc', format='ascii', overwrite=True)
# new_table.write(f'/home/sameer/Thesis_Phase_1/VPfit/PG0003+158_cont_norm.asc', format='ascii', overwrite=True)

# data['CONT_FLUX']=cont_col
# data['NORMALISED_FLUX']=norm_flux_col
# data['NORMALISED_ERROR']=norm_err_col
# data.add_columns([cont_col,norm_flux_col,norm_err_col])
# data.write(file,format='fits', overwrite='True')

quit()

low_cont=cont*0.97
cont_norm_flux=Column(name='FLUX',data=flux/low_cont)
err_low=Column(name='ERROR',data=err/low_cont)

tab=Table()
tab.add_columns([wave,cont_norm_flux,err_low])
tab.write(f'Data/{file[:-5]}_low_cont_norm.fits', format='fits', overwrite=True)
tab.write(f'Data/{file[:-5]}_low_cont_norm.asc', format='ascii', overwrite=True)

hi_cont=cont*1.03
cont_norm_flux=Column(name='FLUX',data=flux/hi_cont)
err_hi=Column(name='ERROR',data=err/hi_cont)

tab=Table()
tab.add_columns([wave,cont_norm_flux,err_hi])
tab.write(f'Data/{file[:-5]}_high_cont_norm.fits', format='fits', overwrite=True)
tab.write(f'Data/{file[:-5]}_high_cont_norm.asc', format='ascii', overwrite=True)
