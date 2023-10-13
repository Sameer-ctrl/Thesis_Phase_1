
# from curses.ascii import isupper
from astropy.io import ascii
from numpy import *
import matplotlib.pyplot as plt
import os

plt.style.use('my_style.mpl')

qso='3c263'
z_abs=0.140723

qso_list=loadtxt('Data/qso_list.txt',dtype=str)
qso_dict=dict(zip(qso_list[:,1],qso_list[:,0]))
qso_label=qso_dict[qso]

param_line=['OVI      0.1407564505    0.0000092985    25.92999      3.94350  13.627117    0.037619',
            'HI1      0.1407019306    0.0000021890    86.81166     10.11454  13.486395    0.064776', 
            'HI2      0.1407546589    0.0000017485    28.38333      0.69949  14.494759    0.022348']

param_line_Shull=['OVI   0 6 33 12 13.60 0.09',
                  'HI1   -7 7 87 15 13.47 0.10', 
                  'HI2    7 1 28 1  14.51 0.03']

param={}
param_Shull={}

keys=['OVI','HI 1','HI 2']

for i,key in enumerate(keys):
    line1=param_line[i]
    splitted1=line1.split()
    param[keys[i]]=[float(x) for x in splitted1[1:]]

    line2=param_line_Shull[i]
    splitted2=line2.split()
    param_Shull[keys[i]]=[float(x) for x in splitted2[1:]]

v_z=lambda z : 3e5*(((1+z)**2-1)/((1+z)**2+1))
del_vz=lambda z,z_err: 4*3e5*((1+z)/(((1+z)**2)+1)**2)*z_err

v_abs=v_z(z_abs)


class line():

    def __init__(self,key,Shull=False):

        if Shull==False:

            z, z_err, self.b, self.b_err, self.logN, self.logN_err = param[key]
            self.v=v_z(z)-v_abs
            self.v_err=del_vz(z,z_err)
        
        else:
            self.v, self.v_err, self.b, self.b_err, self.logN, self.logN_err = param_Shull[key]


line_obj=[line(x) for x in keys]
line_obj_Shull=[line(x,Shull=True) for x in keys]

x=linspace(1,len(keys),len(keys))


v=array([i.v for i in line_obj])
v_err=array([i.v_err for i in line_obj])

v_shull=array([i.v for i in line_obj_Shull])
v_err_shull=array([i.v_err for i in line_obj_Shull])

v=zeros(len(keys))
v_err=zeros(len(keys))
v_shull=zeros(len(keys))
v_err_shull=zeros(len(keys))


b=zeros(len(keys))
b_err=zeros(len(keys))
b_shull=zeros(len(keys))
b_err_shull=zeros(len(keys))

logN=zeros(len(keys))
logN_err=zeros(len(keys))
logN_shull=zeros(len(keys))
logN_err_shull=zeros(len(keys))

for i in range(len(keys)):

    v[i]=line_obj[i].v
    v_err[i]=line_obj[i].v_err
    v_shull[i]=line_obj_Shull[i].v
    v_err_shull[i]=line_obj_Shull[i].v_err

    b[i]=line_obj[i].b
    b_err[i]=line_obj[i].b_err
    b_shull[i]=line_obj_Shull[i].b
    b_err_shull[i]=line_obj_Shull[i].b_err

    logN[i]=line_obj[i].logN
    logN_err[i]=line_obj[i].logN_err
    logN_shull[i]=line_obj_Shull[i].logN
    logN_err_shull[i]=line_obj_Shull[i].logN_err


fig=plt.figure(figsize=(30,6))

plt.subplot(1,3,1)

plt.title('Parameter : v (0 at z = z_abs)')

plt.scatter(x,v,label='Our',color='red')
plt.fill_between(x,v-v_err,v+v_err,alpha=0.2,facecolor='red')

plt.scatter(x,v_shull,label='Shull et al.',color='green')
plt.fill_between(x,v_shull-v_err_shull,v_shull+v_err_shull,alpha=0.2,facecolor='green')

plt.xlabel('line')
plt.ylabel('v (km/s)')
plt.xticks(x,keys)


plt.subplot(1,3,2)

plt.title('Parameter : b')

plt.scatter(x,b,label='Our',color='red')
plt.fill_between(x,b-b_err,b+b_err,alpha=0.2,facecolor='red')

plt.scatter(x,b_shull,label='Shull et al.',color='green')
plt.fill_between(x,b_shull-b_err_shull,b_shull+b_err_shull,alpha=0.2,facecolor='green')

plt.xlabel('line')
plt.ylabel('b (km/s)')
plt.xticks(x,keys)


plt.subplot(1,3,3)

plt.title('Parameter : logN')

plt.scatter(x,logN,label='Our',color='red')
plt.fill_between(x,logN-logN_err,logN+logN_err,alpha=0.2,facecolor='red')

plt.scatter(x,logN_shull,label='Shull et al.',color='green')
plt.fill_between(x,logN_shull-logN_err_shull,logN_shull+logN_err_shull,alpha=0.2,facecolor='green')

plt.xlabel('line')
plt.ylabel(r'logN (cm$^{-2}$)')
plt.xticks(x,keys)

plt.legend(bbox_to_anchor=(0.93,0.95),bbox_transform=plt.gcf().transFigure, loc='center',ncols=1)
plt.suptitle(f'$\mathbf{{{qso_label} \ (z_{{abs}}={z_abs:.6f})}}$')
plt.savefig(f'Files_n_figures/Shull_et_al_comparison_{qso}.png')
# plt.show()



quit()


import subprocess

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

vpfit = subprocess.Popen(['vpfit'], stdin=subprocess.PIPE, std_err=subprocess.PIPE, text=True,shell=True)
output, _errors = vpfit.communicate(input=ip)


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



