from astropy.io import fits,ascii
from matplotlib import table
from numpy import *
import matplotlib.pyplot as plt
from astropy.table import Table
from roman import toRoman
from scipy.interpolate import interp2d,interp1d
from roman import toRoman


plt.style.use('Files_n_figures/my_style.mpl')


hdu=fits.open('Data/component_II_nH_Z_const_T_col_density_param.fits')
data=Table(hdu[1].data)

log_nH=data['log_nH']
log_Z=data['log_Z']
col_den_HII=log10(data['H+'])

nH=[]
Z=[]
col_den=[]

for i in range(len(data)):

    if log_nH[i]>=-3 and 1>=log_Z[i]>=-0.4:
        nH.append(log_nH[i])
        Z.append(log_Z[i])
        col_den.append(col_den_HII[i])


plt.figure()

ax=plt.axes(projection='3d')

ax.scatter(nH,Z,col_den)
# ax.set_xlim(left=-3)
# ax.set_ylim(bottom=-0.4)
# plt.xlim(-3,0)
# plt.ylim(-2,0)
ax.set_xlabel('nH')
ax.set_ylabel('Z')
ax.set_zlabel('HII')

plt.figure()
plt.scatter(nH,col_den,c=Z,cmap='jet')
plt.colorbar()
plt.title('nH')

plt.figure()
plt.scatter(Z,col_den,c=nH,cmap='jet')
plt.colorbar()
plt.title('Z')

plt.show()

quit()

hdu=fits.open('Data/component_II_nH_const_T_col_density_param.fits')
data=Table(hdu[1].data)

nH=data['log_nH']
ions=['Si+', 'Si+2','C+', 'C+2']
ions_label=[f'Si {toRoman(2)}', f'Si {toRoman(3)}',f'C {toRoman(2)}', f'C {toRoman(3)}']

upper_col_den=[12.4,12.3,13.3,'NA']

upper_lim=dict(zip(ions,upper_col_den))

N_OVI=14.26
T_OVI=5.29
log_Zref=-1

col_den_OVI=log10(data['O+5'])

Z=N_OVI-col_den_OVI+log_Zref


ls=['dashed','dashdot','dotted',(0, (3, 5, 1, 5, 1, 5)),'solid']


for i,ion in enumerate(ions):
    plt.plot(nH,log10(data[ion])+Z-log_Zref,ls=ls[i],label=f'{ions_label[i]} ({upper_lim[ion]})')
    plt.scatter(nH,log10(data[ion])+Z-log_Zref)


# plt.scatter(nH,Z)
# plt.plot(nH,Z,ls='--')

# plt.ylabel(r'$\mathbf{log \ [Z]}$',labelpad=15)
plt.ylabel(r'$\mathbf{log \ [N \ {cm}^{-2}]}$',labelpad=15)
plt.xlabel(r'$\mathbf{log \ [n_H \ ({cm}^{-3})]}$',labelpad=15)
plt.legend()
plt.show()


quit()

# mask=T==5.30

# nH_mask=nH[mask]
# col_den_OVI_mask=col_den_OVI[mask]

# plt.scatter(nH_mask,col_den_OVI_mask)
# plt.plot(nH_mask,col_den_OVI_mask,ls='--')
# plt.xlim(right=0.6)
# plt.xlabel(r'$\mathbf{log \ [n_H \ {cm}^{-3}]}$',labelpad=15)
# plt.ylabel(r'$\mathbf{log \ [N \ {cm}^{-2}]}$',labelpad=15)
# plt.show()
# quit()



f=interp1d(T,col_den_OVI,kind='cubic')

N_OVI=14.26

T_plot=linspace(4,7,1000)
interp_col_den=f(T_plot)
log_Zref=0

Z=arange(-1.13,-1.08,0.01)

ax=plt.axes()
axin=ax.inset_axes([0.47,0.41,0.07,0.2])

for z in Z:
    ax.plot(T_plot,interp_col_den+round(z,2),label=f'log Z = {round(z,2)}',ls='--')
    axin.plot(T_plot,interp_col_den+round(z,2),label=f'{round(z,2)}',ls='--')
    

plt.plot(T_plot,interp_col_den,label='Solar metallicity')

axin.vlines(5.29,ymin=1,ymax=N_OVI,lw=3,color='black')
axin.hlines(N_OVI,xmin=3,xmax=5.29,lw=3,color='black')
plt.vlines(5.29,ymin=1,ymax=N_OVI,lw=3,color='black')
plt.hlines(N_OVI,xmin=3,xmax=5.29,lw=3,color='black')
axin.set_xlim(5.2898,5.2902)
axin.set_ylim(14.23,14.28)
axin.set_xticks([])
axin.set_yticks([])
ax.indicate_inset_zoom(axin)
plt.ylim(bottom=10)
plt.legend()
plt.xlim(left=3.8)
plt.xlabel(r'$\mathbf{log \ [T \ (k)]}$',labelpad=15)
plt.ylabel(r'$\mathbf{log \ [N \ {cm}^{-2}]}$',labelpad=15)
plt.show()





'-----------------------'
quit()

# hdu=fits.open('Data/component_II_nH_Z_col_density_param.fits')
# data=Table(hdu[1].data)

# log_nH=data['log_nH']
# log_Z=data['log_Z']

# col_den_OVI=data['O+5']

# Z=-0.30
# N_OVI=14.26

# mask=log_Z==Z

# data=data[mask]

# log_nH=data['log_nH']
# col_den_OVI=data['O+5']

# plt.plot(log_nH,log10(col_den_OVI))
# plt.hlines(N_OVI,xmin=-5,xmax=-2.4145,ls='--',color='black')
# plt.vlines(-2.4145,ymin=14,ymax=N_OVI,ls='--',color='black')
# plt.ylim(bottom=14.1)
# plt.show()


'-----------------------'
# quit()


data=ascii.read('Data/gnat_sternberg_O_CIE.txt')

t=data['log_T']
f_OVI=data['ionf_6']

f_int=interp1d(t,f_OVI,kind='quadratic')

Z_Zsol=0.1
nH=-3.995
N_HI=14.13
N_OVI=14.26
l=10**(N_HI-nH)   #cm


def N_O(x):

    log_nO=x+nH+log10(0.078/91.2)
    col_den=log_nO+N_HI-nH

    return col_den

log_Z=-0.31

ax=plt.axes()
axin=ax.inset_axes([0.3,0.6,0.2,0.2])

for x in arange(4.28,4.32,0.01):
    x=round(x,2)
    ax.scatter(t,log10(f_OVI)+N_O(x)+log_Z-log10(Z_Zsol),label=f'{x}')
    axin.scatter(t,log10(f_OVI)+N_O(x)+log_Z-log10(Z_Zsol),label=f'{x}')

axin.vlines(5.29,ymin=1,ymax=N_OVI,ls='--',color='black')
axin.hlines(N_OVI,xmin=4,xmax=5.29,ls='--',color='black')

plt.vlines(5.29,ymin=1,ymax=N_OVI,ls='--',color='black')
plt.hlines(N_OVI,xmin=4,xmax=5.29,ls='--',color='black')
plt.legend()
plt.xlim(left=4.7)
plt.xlabel(r'$\mathbf{log \ [T \ (k)]}$',labelpad=15)
plt.ylabel(r'$\mathbf{log \ [N \ {cm}^{-2}]}$',labelpad=15)
axin.set_xlim(5.2898,5.2902)
axin.set_ylim(14.23,14.30)
axin.set_xticks([])
axin.set_yticks([])
ax.indicate_inset_zoom(axin)
plt.show()

'-----------------------'


# hdu=fits.open('Data/component_II_nH_col_density_param.fits')
# data=Table(hdu[1].data)

# mask=data['H']!=0
# data=data[mask]

# nH=data['log_nH']
# ions=['Si+', 'Si+2','C+', 'C+2', 'O+5']
# ions_label=[f'Si {toRoman(2)}', f'Si {toRoman(3)}',f'C {toRoman(2)}', f'C {toRoman(3)}',f'O {toRoman(6)}']

# obs_OVI=14.26   

# ls=['dashed','dashdot','dotted',(0, (3, 5, 1, 5, 1, 5)),'solid']

# for i,ion in enumerate(ions):
#     plt.plot(nH,log10(data[ion]),label=ions_label[i],ls=ls[i])

# plt.hlines(obs_OVI,xmin=-6.5,xmax=-4.3729,color='black',ls='--')
# plt.vlines(-4.3729,ymin=1.35,ymax=obs_OVI,color='black',ls='--')
# plt.legend()
# plt.xlim(left=-5.25)
# plt.ylim(bottom=5.22)
# plt.xlabel(r'$log \ [n_H \ ({cm}^{-3})]$',labelpad=15)
# plt.ylabel(r'$log \ [N \ ({cm}^{-2})]$',labelpad=15)
# # plt.show()

# quit()

# bHI=[62.49,2.92]
# bOVI=[29.63,2.04]

# T=(16/(15*(0.129**2)))*(bHI[0]**2-(bOVI[0]**2))
# T_low=(16/(15*(0.129**2)))*((bHI[0]-bHI[1])**2-((bOVI[0]+bOVI[1])**2))
# T_high=(16/(15*(0.129**2)))*((bHI[0]+bHI[1])**2-((bOVI[0]-bOVI[1])**2))

# hdu=fits.open('Data/component_II_nH_T_col_density_param.fits')
# data=Table(hdu[1].data)

# log_nH=data['log_nH']
# log_T=data['log_T']

# nH_ind={}

# for n in linspace(-5.0,0,6):
#     x=[]
#     for i,j in enumerate(log_nH):
#         if j==n:
#             x.append(i)

#     nH_ind[n]=[min(x),max(x)]

# print(nH_ind)

# def nH_ind(nH):

#     x=[]

#     for i,j in enumerate(log_nH):
#         if j==nH:
#             x.append(i)

#     return [min(x),max(x)]

# def interp_func_T(nH):

#     ind=nH_ind(nH)

#     x=[]

#     for j in data['O+6'][ind[0]:ind[1]+1].value:
#         if j==0:
#             x.append(10**(-30))
        
#         else:
#             x.append(j)

#     log_col_den=log10(x)

#     f=interp1d(log_T[ind[0]:ind[1]+1],log_col_den,kind='linear')
# log_T_plot=(linspace(4,7,500))

# f3=interp_func_T(-3)(log_T_plot)
# f2=interp_func_T(-2)(log_T_plot)
# f1=interp_func_T(-1)(log_T_plot)



# plt.plot(log_T_plot,f3,label='-3')
# plt.plot(log_T_plot,f2,label='-2')
# plt.plot(log_T_plot,f1,label='-1')

# plt.legend()
# plt.show()
    # return f

# log_T_plot=(linspace(4,7,500))

# f3=interp_func_T(-3)(log_T_plot)
# f2=interp_func_T(-2)(log_T_plot)
# f1=interp_func_T(-1)(log_T_plot)



# plt.plot(log_T_plot,f3,label='-3')
# plt.plot(log_T_plot,f2,label='-2')
# plt.plot(log_T_plot,f1,label='-1')

# plt.legend()
# plt.show()








