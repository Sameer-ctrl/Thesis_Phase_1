from astropy.io import fits
from numpy import *
import matplotlib.pyplot as plt
from astropy.table import Table
from roman import toRoman
from scipy.interpolate import interp2d,interp1d

plt.style.use('Files_n_figures/my_style.mpl')


hdu=fits.open('Data/component_II_nH_col_density_param.fits')
data=Table(hdu[1].data)

mask=data['H']!=0
data=data[mask]

nH=data['log_nH']
ions=['Si+', 'Si+2','C+', 'C+2', 'O+5']
ions_label=[f'Si {toRoman(2)}', f'Si {toRoman(3)}',f'C {toRoman(2)}', f'C {toRoman(3)}',f'O {toRoman(6)}']

obs_OVI=14.26   

ls=['dashed','dashdot','dotted',(0, (3, 5, 1, 5, 1, 5)),'solid']

for i,ion in enumerate(ions):
    plt.plot(nH,log10(data[ion]),label=ions_label[i],ls=ls[i])

plt.hlines(obs_OVI,xmin=-6.5,xmax=-4.3729,color='black',ls='--')
plt.vlines(-4.3729,ymin=1.35,ymax=obs_OVI,color='black',ls='--')
plt.legend()
plt.xlim(left=-5.25)
plt.ylim(bottom=5.22)
plt.xlabel(r'$log \ [n_H \ ({cm}^{-3})]$',labelpad=15)
plt.ylabel(r'$log \ [N \ ({cm}^{-2})]$',labelpad=15)
plt.show()

quit()

bHI=[62.49,2.92]
bOVI=[29.63,2.04]

T=(16/(15*(0.129**2)))*(bHI[0]**2-(bOVI[0]**2))
T_low=(16/(15*(0.129**2)))*((bHI[0]-bHI[1])**2-((bOVI[0]+bOVI[1])**2))
T_high=(16/(15*(0.129**2)))*((bHI[0]+bHI[1])**2-((bOVI[0]-bOVI[1])**2))

hdu=fits.open('Data/component_II_nH_T_col_density_param.fits')
data=Table(hdu[1].data)

log_nH=data['log_nH']
log_T=data['log_T']

# nH_ind={}

# for n in linspace(-5.0,0,6):
#     x=[]
#     for i,j in enumerate(log_nH):
#         if j==n:
#             x.append(i)

#     nH_ind[n]=[min(x),max(x)]

# print(nH_ind)

def nH_ind(nH):

    x=[]

    for i,j in enumerate(log_nH):
        if j==nH:
            x.append(i)

    return [min(x),max(x)]

def interp_func_T(nH):

    ind=nH_ind(nH)

    x=[]

    for j in data['O+6'][ind[0]:ind[1]+1].value:
        if j==0:
            x.append(10**(-30))
        
        else:
            x.append(j)

    log_col_den=log10(x)

    f=interp1d(log_T[ind[0]:ind[1]+1],log_col_den,kind='linear')

    return f

log_T_plot=(linspace(4,7,500))

f3=interp_func_T(-3)(log_T_plot)
f2=interp_func_T(-2)(log_T_plot)
f1=interp_func_T(-1)(log_T_plot)



plt.plot(log_T_plot,f3,label='-3')
plt.plot(log_T_plot,f2,label='-2')
plt.plot(log_T_plot,f1,label='-1')

plt.legend()
plt.show()





quit()

def interp_func():

    x=[]

    for j in data['O+6'].value:
        if j==0:
            x.append(10**(-30))
        
        else:
            x.append(j)

    log_col_den=log10(x)
    f=interp2d(log_nH,log_T,log_col_den,kind='cubic')

    return f


interp_func_OVI=interp_func()

col_den_med=interp_func_OVI(log_nH,log10(T))
col_den_low=interp_func_OVI(log_nH,log10(T_low))
col_den_high=interp_func_OVI(log_nH,log10(T_high))


plt.figure()

plt.plot(log_nH,col_den_med,label='med',color='red')
plt.plot(log_nH,col_den_low,label='low',color='green',ls='--')
plt.plot(log_nH,col_den_high,label='high',color='blue',ls='--')

plt.legend()
plt.show()