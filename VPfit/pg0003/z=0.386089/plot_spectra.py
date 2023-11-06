from astropy.io import ascii
import matplotlib.pyplot as plt


qso=input('Enter QSO name \n')
data=ascii.read(f'{qso}_cont_norm.asc')

wave=data['WAVE']
flux=data['FLUX']
err=data['ERROR']

plt.step(wave,flux)
plt.step(wave,err)
plt.ylim(-0.3,2.5)
plt.show()