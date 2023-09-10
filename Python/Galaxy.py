from astropy.io import fits
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from numpy import *
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord

plt.style.use('my_style.mpl')

cosmo=FlatLambdaCDM(H0=69.6,Om0=0.3,Tcmb0=2.725)
z_abs=0.347


def L_by_Lstar(m_I,M_star):

    dL=cosmo.luminosity_distance(z_abs)  # Mpc
    M_star_I=5*log10(cosmo.h)+M_star
    M_I=m_I-(5*log10(dL.value*1e6))+5
    L_by_Lstar_val=10**((M_star_I-M_I)/2.5)

    print(L_by_Lstar_val)

L_by_Lstar(17.8,-21.64)
# quit()

v_abs=3e5*(((1+z_abs)**2-1)/((1+z_abs)**2+1))
z1=sqrt((1+((v_abs+1000)/3e5))/(1-((v_abs+1000)/3e5)))-1
z2=sqrt((1+((v_abs-1000)/3e5))/(1-((v_abs-1000)/3e5)))-1

dA_scale=cosmo.kpc_proper_per_arcmin(z_abs)
# print(dA_scale*60)

data=loadtxt('Data/gal_data.txt',dtype=str)


ra=data[:,0].astype(float)
dec=data[:,1].astype(float)
z=data[:,2].astype(float)
v=3e5*(((1+z)**2-1)/((1+z)**2+1))-v_abs


ra_qso,dec_qso=[1.4968167,16.163614]

abs_sys_coord=SkyCoord(ra_qso*u.degree,dec_qso*u.degree)
gal_coords=SkyCoord(ra*u.degree,dec*u.degree)


def circle(r,x0,y0,label):

    theta=linspace(0,2*pi,1000)
    x=x0+(r*cos(theta))
    y=y0+(r*sin(theta))

    plt.plot(x,y,ls='--',label=label)

print(abs(v))

plt.figure(figsize=(16,13))

plt.scatter((ra-ra_qso)*60,(dec-dec_qso)*60,s=50,c=abs(v),cmap='jet')

# for i in range(len(ra)):
#     plt.text((ra[i]-ra_qso)*60,(dec[i]-dec_qso)*60,s=f'{i+1}')

cb=plt.colorbar(label=r'$\mathbf{\Delta V \ (km \ s^{-1})}$')
cb.ax.tick_params(labelsize=25)
cb.set_label(r'$\mathbf{\Delta V \ (km \ s^{-1})}$',labelpad=40)
circle(500/dA_scale.value,0,0,r'$\mathbf{500 \ kpc}$')
circle(1000/dA_scale.value,0,0,r'$\mathbf{1 \ Mpc}$')
# circle(5000/dA_scale.value,0,0,'5 Mpc')
plt.scatter(0,0,s=200,c='black',marker='X')
plt.legend(fontsize=30)
plt.xlim(-15,15)
plt.ylim(-15,15)
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
plt.xlabel(r'$\mathbf{\alpha-{\alpha}_{QSO} \ [arc \ min]}$',fontsize=30,labelpad=15)
plt.ylabel(r'$\mathbf{\delta-{\delta}_{QSO} \ [arc \ min]}$',fontsize=30)
# plt.savefig('galaxy_environment.png')
plt.show()


# plt.colorbar(label='separation')
# plt.scatter(0,0,s=200,c='black')

# plt.show()











'Luminosity limit'

# dL=cosmo.luminosity_distance(z_abs)  # Mpc

# Mstar_AB_r=5*log10(cosmo.h)-21.64    # AB magnitude
# Mstar_r=Mstar_AB_r-0.055             #  Johnson
# m_r=17.7
# M_r=m_r-(5*log10(dL.value*1e6))+5

# L_by_Lstar=10**((Mstar_r-M_r)/2.5)

# print(L_by_Lstar)

