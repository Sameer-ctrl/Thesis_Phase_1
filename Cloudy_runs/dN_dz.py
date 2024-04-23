from  numpy import *
from scipy.constants import G, m_e,m_p,speed_of_light,parsec



rest_wave=1215.6701
wave=array([1230,1556.5])
z=(wave-rest_wave)/rest_wave

# print(z[1]-z[0])

z_qso=0.297

v_z=lambda z : 3e5*(((1+z)**2-1)/((1+z)**2+1))  # v at z
z_v=lambda v : sqrt((1+((v)/3e5))/(1-((v)/3e5)))-1      # z at v

v_qso=v_z(z_qso)
v_z1=v_z(z[1])

# print(v_qso-v_z1)

z_lim=z_v(v_qso-5000)
 

# print((z_lim+1)*rest_wave)


wave=array([1230,1556.5])
z=(wave-rest_wave)/rest_wave
dz_unbl=0.238
NH=array([18.49,18.23,18.58,18.51,18.80,18.70])

# print(log10(sum(10**NH)))

zmax=z[1]
zmin=z[0]
dzb=zmax-zmin-dz_unbl

dX=0.5*((1+zmax-(dzb/2))**2-(1+zmin+(dzb/2))**2)    # ref Richter et al. 2006 
print(dX)


mu=1.3
mH=m_p+m_e
Ho=70
rho_c=(3*Ho**2)/(8*pi*G)

Omega_b_BLA1=((mu*mH*Ho)/(rho_c*speed_of_light))*sum(10**NH)*(1/dX)

print(Omega_b_BLA1*(1e7*parsec)*100)
