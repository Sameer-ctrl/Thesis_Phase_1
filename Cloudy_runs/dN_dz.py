from  numpy import *
from scipy.constants import G, m_e,m_p,speed_of_light,parsec
import os
from astropy.io import ascii,fits
import matplotlib.pyplot as plt
from astropy.table import Table


# rest_wave=1215.6701
# wave=array([1230,1556.5])
# z=(wave-rest_wave)/rest_wave

# # print(z[1]-z[0])

# z_qso=0.297

# v_z=lambda z : 3e5*(((1+z)**2-1)/((1+z)**2+1))  # v at z
# z_v=lambda v : sqrt((1+((v)/3e5))/(1-((v)/3e5)))-1      # z at v

# v_qso=v_z(z_qso)
# v_z1=v_z(z[1])

# # print(v_qso-v_z1)

# z_lim=z_v(v_qso-5000)
 

# # print((z_lim+1)*rest_wave)


# wave=array([1230,1556.5])
# z=(wave-rest_wave)/rest_wave
# dz_unbl=0.238
# NH=array([18.49,18.23,18.58,18.51,18.80,18.70])

# # print(log10(sum(10**NH)))

# zmax=z[1]
# zmin=z[0]
# dzb=zmax-zmin-dz_unbl

# dX=0.5*((1+zmax-(dzb/2))**2-(1+zmin+(dzb/2))**2)    # ref Richter et al. 2006 
# print(dX)


# mu=1.3
# mH=m_p+m_e
# Ho=70
# rho_c=(3*Ho**2)/(8*pi*G)

# Omega_b_BLA1=((mu*mH*Ho)/(rho_c*speed_of_light))*sum(10**NH)*(1/dX)

# print(Omega_b_BLA1*(1e7*parsec)*100)



qso='h1821'

file=f'../Python/Data/IGM_Danforth_Data/Cont_norm_spectra/{qso}_cont_norm.asc'

c=3e5

data=ascii.read(file)
wave=data['WAVE']
flux=data['FLUX']

data_linelist=ascii.read(f'../Python/Data/IGM_Danforth_Data/Linelist/{qso}_linelist.txt')
col_names_linelist=['WAVELENGTH','LINE_ID','z_ABS','SIGLEVEL','S/N','EQW','EQW_ERR','BVALUE','BVALUE_ERR','LOGN_lower_limit','LOGN','LOGN_ERR','FLAG_FIT','LCOD','FLAG_ISM','FLAG_AGN']
data_linelist.rename_columns(data_linelist.colnames,col_names_linelist)

wave_abs=data_linelist['WAVELENGTH']
eqw_abs=data_linelist['EQW']
b_abs=data_linelist['BVALUE']
z_abs=data_linelist['z_ABS']

line_id=data_linelist['LINE_ID']

data_rest=loadtxt('../Python/Data/rest_wave.txt',dtype=str)

line_atom=data_rest[:,0]
wave_rest=data_rest[:,1].astype(float)

rest_wave_dict={}

for i in range(len(line_atom)):
    rest_wave_dict.update({line_atom[i]:wave_rest[i]})




def rest_wave(lines):

    r_wave=zeros(len(lines))
    

    for i,line in enumerate(lines):

        if z_abs[i]!=0:
        
            split=line.split()

            if split[0][:2]=='Ly':
                key='_'.join(['HI',split[1]])

            else:
                key='_'.join(split)


            if key!='NULL':
                r_wave[i]=rest_wave_dict[key]

            else:
                r_wave[i]=rest_wave_dict['HI_1215']

        else:
            r_wave[i]=rest_wave_dict['HI_1215']

    return r_wave

rest_wavelength=rest_wave(line_id)
print(2*sqrt(log(2)))

plt.step(wave,flux)
plt.errorbar(wave_abs,ones(len(wave_abs)),yerr=0,xerr=2*sqrt(log(2))*(b_abs/c)*rest_wavelength,fmt='o',capsize=3)
plt.ylim(-0.2,2)
plt.show()