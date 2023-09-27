import os
from astropy.io import ascii
from astropy.table import Table
from numpy import *


b='Spectra'

files=os.listdir(f'Data/IGM_Danforth_Data/{b}')

for file in files:
    ind=[]
    for i,a in enumerate(file):
        if a=='_':
            ind.append(i)
        
    qso=file[ind[3]+1:ind[4]]
    os.rename(f'Data/IGM_Danforth_Data/{b}/{file}',f'Data/IGM_Danforth_Data/{b}/{qso}_spec.fits')
