from curses.ascii import isupper
from astropy.io import ascii
from astropy.io import fits
from astropy.table import Table
from numpy import *
import subprocess
import matplotlib.pyplot as plt
import os


'finding contamination missed by Danforth'


qso='pg1216'
# z_abs=0.08092
# cen_wave_rest=1550.7812

# x=(z_abs+1)*cen_wave_rest
x=1224.3

data=loadtxt('../Python/Data/rest_wave.txt',dtype=str)

data_sys=ascii.read(f'../Python/Data/IGM_Danforth_Data/Systems/{qso}_igm-systems.txt')

z_sys=unique(data_sys['col1'])
# print(z_sys)
line_atom=data[:,0]
wave_atom=data[:,1].astype(float)

for i,w in enumerate(wave_atom):
    z=(x-w)/w

    arg_min=argmin(abs(z_sys-z))
    close_z=z_sys[arg_min]

    if 0<=z<=0.9:

        if abs(z-close_z)*1e6<1000:

            print(f'{line_atom[i]} : {z:.6f}  : {close_z:.6f}  : {abs(z-close_z)*1e6:.6f}')


quit()
qso='pg0832'

file=f'Data/IGM_Danforth_Data/Spectra/{qso}_spec.fits'


hdu=fits.open(file)
data=Table(hdu[1].data)
wave=data['WAVE']
flux=data['FLUX']

# z=0.077419

# cen_wave_obs=1393.76018*(1+z)

plt.step(wave,flux,label='SiIV 1393')
plt.step(wave*(1393.76018/1402.77291),flux,label='SiIV 1402',ls='--')
plt.legend()
plt.show()






quit()

# def tick_pos():

#     with open('CIV_1548.txt') as f:
#         data=f.readlines()

#     i=0 
#     tick_wave=[]

#     for line in data:
        
#         if line[0]=='!':
#             i+=1

#             if i>2:
#                 tick_wave.append(float(line[:-1].split()[1]))

#     return array(tick_wave)


# plt.vlines(tick_pos(),[0],[1])
# plt.show()

# quit()

# files=os.listdir('pg1216/HI_chunks')

os.chdir('h1821/z=0.224981')

rest_wave=[1025,1215,972,949,937,930,926,923,920,919,918,917,916,'Ly14','Ly15','Ly16','Ly17']
comp=3

for i,wave in enumerate(rest_wave):

    if i<9:
        file=f'vpfit_chunk00{i+1}.txt'

    else:
        file=f'vpfit_chunk0{i+1}.txt'

    if comp==0:
    
        if i<13:
            os.rename(file,f'HI_{wave}.txt')
            # os.rename(file,f'HI_{wave}_{comp}.txt')
        
        else:
            os.rename(file,f'{wave}.txt')
            # os.rename(file,f'{wave}_{comp}.txt')
    
    else:

        if i<13:
            os.rename(file,f'HI_{wave}_{comp}.txt')
        
        else:
            os.rename(file,f'{wave}_{comp}.txt')








quit()
# a='HI_1215'
a='HeeeII_1206'

def line_label(line):

    for i,s in enumerate(line):
        if s.isupper():
            if i==0:
                pass
            else:
                break

    ion=line[:i]
    n,wave=line[i:].split('_')

    return ion,n,wave


quit()

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
