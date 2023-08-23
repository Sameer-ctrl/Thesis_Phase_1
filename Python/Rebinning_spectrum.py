from numpy import *
from astropy.io import fits,ascii
from astropy.table import Table, Column
import astropy.units as u
import matplotlib.pyplot as plt
from linetools.spectra.xspectrum1d import XSpectrum1D
from scipy.integrate import simpson


plt.style.use('my_style.mpl')

hdu=fits.open('Data/PG0003+158.fits')
data=Table(hdu[1].data)
# data_a=ascii.read('Data/spec_PG0003+158_v3.dat')


wave=data['WAVE'][data['WAVE']>=1132.7]
flux=data['FLUX'][data['WAVE']>=1132.7]
err=data['ERROR'][data['WAVE']>=1132.7]

# wave_a=data_a['WAVE']
# flux_a=data_a['FLUX']
# err_a=data_a['ERROR']
# cont_a=data_a['CONT']


spec=XSpectrum1D.from_tuple((wave, flux, err), verbose=False)
rebin_wave=arange(wave[0],wave[-1],0.03669)[:-1]*u.AA
rebin_spec=spec.rebin(rebin_wave,do_sig=True,grow_bad_sig=True)
rebin_flux=rebin_spec.flux
rebin_err=rebin_spec.sig

tab=Table()

wave_col=Column(name='WAVE',data=rebin_wave)
flux_col=Column(name='FLUX',data=rebin_flux)
err_col=Column(name='ERROR',data=rebin_err)

tab.add_columns([wave_col,flux_col,err_col])
tab.write(f'Data/PG0003+158_rebinned.fits', format='fits', overwrite=True)
tab.write(f'Data/PG0003+158_rebinned.asc', format='ascii', overwrite=True)




# plt.step(wave,err,label='original')
# plt.step(rebin_wave,rebin_spec.flux,label='rebinned')
# plt.step(rebin_wave,rebin_spec.sig,label='rebinned')
# plt.step(wave_a,err_a,label='anand')
# plt.legend()
# plt.show()
