from  numpy import *
from scipy.constants import *
import os
from astropy.io import ascii,fits
import matplotlib.pyplot as plt
from astropy.table import Table,Column

mH=(proton_mass+electron_mass)

v_z=lambda z : 3e5*(((1+z)**2-1)/((1+z)**2+1))  # v at z
z_v=lambda v : sqrt((1+((v)/3e5))/(1-((v)/3e5)))-1      # z at v

def redshift_path(qso,wave_min=1220,v_lim=5000):

    file_systems=open(f'Data/IGM_Danforth_Data/Systems/{qso}_igm-systems.txt','r')
    z_em=float(file_systems.readlines()[16].split(' ')[1])

    data=ascii.read(f'Data/IGM_Danforth_Data/Cont_norm_spectra/{qso}_cont_norm.asc')
    wave=data['WAVE']

    excluded_wave=ascii.read(f'Data/IGM_Danforth_Data/Excluded_wavelengths/{qso}_excluded_wavelength.asc')

    wave_l=excluded_wave['WAVE1']
    wave_r=excluded_wave['WAVE2']

    z_lim=z_v(v_z(z_em)-v_lim)

    rest_wave=1215.6701
    wave_max=min([(1+z_lim)*rest_wave,wave[-1]])

    dzb=0

    for i in range(len(excluded_wave)):

        if wave_l[i] >= wave_min:
            dzb+=(wave_r[i]-wave_l[i])/rest_wave
        
        elif wave_r[i] >= wave_min > wave_l[i]:
            dzb+=(wave_r[i]-wave_min)/rest_wave

    dzb=round(dzb,3)

    zmax=round((wave_max-rest_wave)/rest_wave,3)
    zmin=round((wave_min-rest_wave)/rest_wave,3)

    dz=zmax-zmin
    dz_unblocked=round(dz-dzb,3)

    delta_X=round(0.5*(((1+zmax-(dzb/2)))**2-((1+zmin+(dzb/2)))**2),3)
    # print(delta_X,dz_unblocked,dzb,dz)

    return delta_X,dz_unblocked,dzb,dz

def f_H(b):

    'f_H=(HI+HII)/HI ~ HII/HI'

    logT=log10(mH*((b*1000)**2)*(1/(2*Boltzmann)))
    log_fH=(5.4*logT)-(0.33*(logT**2))-13.9   # ref. Sutherland & Dopita 1993, Richter et al. 2004, Tracing baryons in WHIM using BLAs

    return log_fH

h=0.7
H0=100*h    #km/s/ Mpc
mu=1.3

rho_c=2*(H0**2)*(1/(8*pi*gravitational_constant))

A=(mu*mH*H0)/(rho_c*speed_of_light)*(parsec*1e7)   # cm^2

def omega_BLA(qso,b,N):

    if not isinstance(qso, (list, tuple, type(array([])))):
        qso=array([qso])

    if not isinstance(b, (list, tuple, type(array([])))):

        b=array([b])
        N=array([N])

    if isinstance(b, (list, tuple)):
        b=array(b)
        
    if isinstance(N, (list, tuple)):
        N=array(N)
            
    print(type(b),type(N))
    delta_X=sum([redshift_path(q)[0] for q in unique(qso)])
    fH=f_H(b)
    NH=sum(10**(fH+N))

    return A*(NH/delta_X)


qso=['3c263', 'pks0637', 'pks0637', 'pg1424', 'pg0003', 'pg0003', 'pg0003', 'pg1216', 's135712', '1es1553', 'sbs1108', 'pg1222', 'pg1116', 'h1821', 'h1821', 'pg1121', 'pks0405']
b=array([87.0, 162.0, 46.0, 29.0, 63.0, 40.0, 64.0, 52.0, 46.0, 51.0, 16.0, 52.0, 71.0, 63.0, 84.0, 60.0, 26.0])
N=array([13.49, 13.6, 14.61, 15.44, 14.2, 14.1, 14.17, 15.1, 15.01, 13.88, 15.79, 14.34, 13.6, 13.68, 13.64, 14.34, 13.46])


omega_BLA_all=omega_BLA(qso,b,N)*100
print(omega_BLA_all)
print(mH*((40*1000)**2)*(1/(2*Boltzmann)))

omega_BLA_los=[]

for i in range(len(qso)):
    
    los=[qso[i]]
    b_los=array([b[i]])
    N_los=array([N[i]])

    omega_BLA_los.append(omega_BLA(los,b_los,N_los)*100)

plt.figure()

plt.hist(omega_BLA_los,bins='auto')
plt.vlines(omega_BLA_all,0,7,ls='--',color='red')

plt.figure()
plt.scatter(qso,omega_BLA_los)
plt.hlines(omega_BLA_all,qso[0],qso[-1])

plt.figure()

plt.plot(linspace(0,150,1000),f_H(linspace(0,150,1000)))

plt.show()

qso='pg0003'
