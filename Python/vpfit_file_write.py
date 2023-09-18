from operator import indexOf
from astropy.io import ascii
import os

from numpy import argmin, mean

data=ascii.read('Data/abs_system.csv')
lsf_files=os.listdir('Data/COS_LSF')


def lsf_file_search(wave_range):

    file_list=[]
    wave_diff1=[]
    wave_diff2=[]

    for file in lsf_files:
        wave_lsf=int(file[14:18])

        if wave_range[0]<=wave_lsf<=wave_range[1]:
            file_list.append(file)

        else:
            wave_diff1.append(abs(wave_range[0]-wave_lsf))
            wave_diff2.append(abs(wave_range[1]-wave_lsf))

    if len(file_list)>0:
        
        wave_diff_cen=[]

        for file in file_list:
            wave=int(file[14:18])
            wave_diff_cen.append(abs(wave-mean(wave_range)))
        
        i=argmin(wave_diff_cen)

        return file_list[i]
        
    else:
    
        if min(wave_diff1) > min(wave_diff2):
            return lsf_files[argmin(wave_diff2)]
        
        else:
            return lsf_files[argmin(wave_diff1)]



