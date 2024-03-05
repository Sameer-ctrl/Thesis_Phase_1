from astropy.io import ascii
from astropy.table import Table
import matplotlib.pyplot as plt
from numpy import *
from scipy.integrate import simpson

qso='3c263'
z_abs=0.140756

data=ascii.read(f'../Python/Data/IGM_Danforth_Data/Cont_norm_spectra/{qso}_cont_norm.asc')

flux=data['FLUX']
wave=data['WAVE']
err=data['ERROR']


atom_data=loadtxt('../Python/Data/lines_system_plots.txt',dtype=str)

line_atom=atom_data[:,0]
wave_rest=atom_data[:,1].astype(float)
f_osc=atom_data[:,2].astype(float)

rest_wave={}
f={}

for i in range(len(line_atom)):
    rest_wave.update({line_atom[i]:wave_rest[i]})
    f.update({line_atom[i]:f_osc[i]})

lines=['SiII_1193','SiII_1190']
vlim=[-100,100]


for line in lines:

    f_line=f[line]
    cen_wave_rest=rest_wave[line]
    cen_wave_obs=cen_wave_rest*(1+z_abs)

    v=3*(10**5)*((wave**2-(cen_wave_obs**2))/(wave**2+(cen_wave_obs**2)))

    mask=logical_and(v>vlim[0],v<vlim[1])

    v=v[mask]
    tau_v=-log(flux[mask])

    N=simpson(tau_v,v)/((f_line*cen_wave_rest*2.654E-15))

    print(log10(abs(N)))

# plt.subplot(121)
# plt.step(v_mask,N_v)

# plt.subplot(122)
# plt.step(v,flux)
# plt.xlim(-100,100)
# plt.ylim(-0.1,1.25)
# plt.show()
