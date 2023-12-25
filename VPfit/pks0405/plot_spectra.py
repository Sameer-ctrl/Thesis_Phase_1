from astropy.io import ascii
import matplotlib.pyplot as plt


qso=input('Enter QSO name \n')
data=ascii.read(f'{qso}_cont_norm.asc')
data1=ascii.read(f'{qso}_cont_norm_previous.asc')


wave=data['WAVE']
flux=data['FLUX']
err=data['ERROR']

wave1=data1['WAVE']
flux1=data1['FLUX']
err1=data1['ERROR']

plt.step(wave,flux,label='new')
plt.step(wave1,flux1,label='previous')
# plt.step(wave,err)
plt.legend()
plt.show()