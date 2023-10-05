from numpy import *
from astropy.io import fits,ascii
from astropy.table import Table, Column
import astropy.units as u
import matplotlib.pyplot as plt
from linetools.spectra.xspectrum1d import XSpectrum1D
from scipy.integrate import simpson
import os


plt.style.use('my_style.mpl')

qso='pg1116'

hdu=fits.open(f'Data/IGM_Danforth_Data/Spectra/{qso}_spec.fits')
data=Table(hdu[1].data)

wave=data['WAVE']
flux=data['FLUX']
err=data['ERR']

# plt.step(wave,flux)
# plt.step(wave,err)
# plt.show()
# quit()

'Rebinning'

spec=XSpectrum1D.from_tuple((wave, flux, err), verbose=False)
rebin_wave=arange(wave[0],wave[-1],0.03669)[:-1]*u.AA
rebin_spec=spec.rebin(rebin_wave,do_sig=True,grow_bad_sig=True)
rebin_flux=rebin_spec.flux
rebin_err=rebin_spec.sig


'Continuum fitting'


spec = XSpectrum1D.from_tuple((rebin_wave,rebin_flux,rebin_err))
spec.fit_continuum(kind='QSO', redshift=z)

cont=spec.co

tab=Table()

wave_col=Column(name='WAVE',data=rebin_wave)
flux_col=Column(name='FLUX',data=rebin_flux/cont)
err_col=Column(name='ERROR',data=rebin_err/cont)

tab.add_columns([wave_col,flux_col,err_col])
tab.write(f'Data/IGM_Danforth_Data/Processed_spectra/{qso}_rebinned_cont_norm.asc', format='ascii', overwrite=True)
tab.write(f'../VPfit/{qso}/{qso}_rebinned_cont_norm.asc', format='ascii', overwrite=True)


plt.step(wave,flux,label='original')
plt.step(wave,err,label='error')
plt.plot(rebin_wave,cont,label='continuum')
plt.step(rebin_wave,rebin_flux,label='rebinned')
plt.step(rebin_wave,rebin_err,label='rebinned_error')
plt.legend()
plt.show()
