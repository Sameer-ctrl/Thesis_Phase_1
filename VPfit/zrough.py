from astropy.io import ascii
from numpy import *
import subprocess
import matplotlib.pyplot as plt
import os



n=2

file='fit_HI.asc'
ip=''

if n==1:
    ip=f'f\n\n{file}\n\nas\n\n\n'

else:
    for i in range(1,n+1):
        if i==1:
            ip+=f'f\n\n{file}\n\n\nas\n\n\n\n'
        
        else:
            ip+=f'\nas\nas\n\n\n\n'

vpfit = subprocess.Popen(['vpfit'], stdin=subprocess.PIPE, stderr=subprocess.PIPE, text=True,shell=True)
output, errors = vpfit.communicate(input=ip)


for i in range(n):

    data=loadtxt(f'vpfit_chunk00{i+1}.txt',comments='!')
    wave=data[:,0]
    flux=data[:,1]
    cont=data[:,3]

    plt.subplot(int(ceil(n/3)),3,i+1)
    plt.plot(wave,cont,c='red')
    plt.step(wave,flux,c='green')
    plt.title(f'{i+1}',fontsize='20')


os.system('rm vpfit_chunk*')
plt.show()


# -----------------------------

# b='Spectra'

# files=os.listdir(f'Data/IGM_Danforth_Data/{b}')

# for file in files:
#     ind=[]
#     for i,a in enumerate(file):
#         if a=='_':
#             ind.append(i)
        
#     qso=file[ind[3]+1:ind[4]]
#     os.rename(f'Data/IGM_Danforth_Data/{b}/{file}',f'Data/IGM_Danforth_Data/{b}/{qso}_spec.fits')
