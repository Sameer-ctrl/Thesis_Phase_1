from astropy.io import ascii
import matplotlib.pyplot as plt
import os


files=os.listdir()

for f in files:
    if f[-4:]=='.asc':
        spectra=f
        break

data=ascii.read(spectra)

wave=data['WAVE']
flux=data['FLUX']
err=data['ERROR']

plt.step(wave,flux)
plt.step(wave,err)
plt.ylim(-0.3,2.5)
plt.show()
