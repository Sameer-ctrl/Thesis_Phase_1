from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
import matplotlib.pyplot as plt


hdu=fits.open('Data/PG0003+158_cont_norm.fits')
data=Table(hdu[1].data)

abs_sys=ascii.read('Data/abs_system.csv')

wave=data['WAVE']
flux=data['FLUX']
err=data['ERROR']

# plt.step(wave,flux)
# plt.plot(wave,data['cont_flux'])
# plt.xlim(1615,1660)
# plt.show()

# plt.figure()
# plt.title('Complete Spectrum')
# plt.step(wave,flux*(10**14))
# plt.step(wave,err*(10**14))
# plt.ylim(0,1)
# # plt.xlim(1615,1660)
# plt.xlim(left=1110)

lines=['Lya 1215','Lyb 1025','Lyg 972','OVI 1032','OVI 1038','CIII 977','CII 1036','SiIII 1206','SiII 1260']

# z_abs=0.347576
z_abs=0.347961

wave_abs=abs_sys['WAVELENGTH'][abs_sys['Z_SYS']==z_abs]
abs_line=abs_sys['LINE_ID'][abs_sys['Z_SYS']==z_abs]

int=1.5

def spec_slice(cen_wave):

    wave_sort=[]
    flux_sort=[]
    vr=[]

    for i in range(len(wave)):
        if cen_wave-int < wave[i] < cen_wave+int:
            wave_sort.append(wave[i])
            flux_sort.append(flux[i])
            vr.append((3*10**5)*((wave[i]**2-(cen_wave**2))/(wave[i]**2+(cen_wave**2))))
    
    plt.step(wave_sort,flux_sort)
    # plt.xlim(-300,300)

for i in range(len(wave_abs)):

    plt.subplot(6,4,i+1)
    spec_slice(wave_abs[i],)
    plt.title(f'{abs_line[i]} ({wave_abs[i]})')

plt.suptitle(f'z_abs={z_abs}',fontsize=15)
plt.show()




