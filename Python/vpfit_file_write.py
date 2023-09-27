from operator import indexOf
from astropy.io import ascii
import os
from numpy import argmin, linspace, logical_and, mean, sqrt, hstack,vstack, where, zeros

data=ascii.read('Data/IGM_systems/hlsp_igm_hst_cos_pg0003_g130m-g160m_v3_igm-systems.txt')
lsf_files=os.listdir('Data/COS_LSF')

lines=data['LINE_ID']


def lsf_file_search(lambda1,lambda2):

    cen_wave=mean([lambda1,lambda2])
    wave_diff=[abs(cen_wave-int(file[14:18])) for file in lsf_files]
    i=argmin(wave_diff)

    return lsf_files[i]

z_sys=0.347
v_sep_lim=200

v_sys=3e5*(((1+z_sys)**2-1)/((1+z_sys)**2+1))

z1=sqrt((1+((v_sys-v_sep_lim)/3e5))/(1-((v_sys-v_sep_lim)/3e5)))-1
z2=sqrt((1+((v_sys+v_sep_lim)/3e5))/(1-((v_sys+v_sep_lim)/3e5)))-1

mask=logical_and(data['Z_SYS']>=z1, data['Z_SYS']<=z2)
data=data[mask]
lines_z=data['LINE_ID'] 
z_lines=data['z_ABS']
b=data['BVALUE']
logN=data['LOGN']
wave=data['WAVELENGTH']

ions=[line.split(' ')[0] for line in lines_z]

ion_to_fit='CII'
ind=[]

for i,ion in enumerate(ions):
    if ion==ion_to_fit:
        ind.append(i)

guess_z=z_lines[ind]
guess_b=b[ind]
guess_logN=logN[ind]

init_guess=vstack((guess_z,zeros(len(ind)),guess_b,zeros(len(ind)),guess_logN,zeros(len(ind)))).transpose()

print(init_guess)






