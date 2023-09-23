# %%
from astropy.io import fits,ascii
from matplotlib import table
from numpy import *
import matplotlib.pyplot as plt
from astropy.table import Table
from roman import toRoman
from scipy.interpolate import interp2d,interp1d
from scipy.optimize import fsolve
from roman import toRoman

# %%
plt.style.use('Files_n_figures/my_style.mpl')

hdu=fits.open('Data/component_II_nH_Z_const_T_col_density_param.fits')
data=Table(hdu[1].data)

log_nH=data['log_nH']
log_Z=data['log_Z']
col_den_OVI=log10(data['O+5'])

def ion_label(ion,state):

    return f'{{\\fontsize{{25pt}}{{3em}}\selectfont{{}}$\mathbf{{{ion}}}$}} {{\\fontsize{{17pt}}{{3em}}\selectfont{{}}$\mathbf{{{state}}}$}}'


N_OVI=14.26

ions=['Si+', 'Si+2','C+', 'C+2']
ions_label=[ion_label('Si','III'),ion_label('Si','II'),ion_label('C','III'),ion_label('C','II'),ion_label('O','VI')]

upper_col_den=[12.4,12.3,13.3,'NA']
upper_lim=dict(zip(ions,upper_col_den))

# %%

# for z in arange(-3,-2,0.05):
#     z=round(z,2)
#     mask=log_Z==-z
#     nH=log_nH[mask]

#     plt.plot(nH,log10(data['H+'])[mask],label=f'{z}')

# plt.legend()
# plt.show()


# %%
def Z_solve(nH):

    mask=log_nH==nH

    log_Z_mask=log_Z[mask]
    col_den_OVI_mask=col_den_OVI[mask]

    f=interp1d(log_Z_mask,col_den_OVI_mask,kind='cubic')

    def func(x):
        return f(x)-N_OVI

    Z=fsolve(func,-1)
    col_den_ions_dict={}

    for i in ions:
        col_den_ion=log10(data[i])[mask]
        f_interp=interp1d(log_Z_mask,col_den_ion,kind='cubic')
        col_den_ions_dict[i]=f_interp(Z)

    return Z, col_den_ions_dict

    
nH=arange(-5,0.1,0.1)
Z=[]
# col_den_ions_dict={'Si+':[],'Si+2':[],'C+':[],'C+2':[]}

col_den_ions_dict={'Si+':[],'Si+2':[],'C+':[],'C+2':[]}
# ions_label=[f'Si {toRoman(2)}', f'Si {toRoman(3)}',f'C {toRoman(2)}', f'C {toRoman(3)}']
ls=['dashed','dashdot','dotted',(0, (3, 5, 1, 5, 1, 5)),'solid']


for n in nH:
    n=round(n,2)
    z,col_den_ions=Z_solve(n)
    if n==-3.0:
        print(z)
    Z.append(z)

    for i in ions:
        col_den_ions_dict[i].append(col_den_ions[i])


plt.figure(figsize=(16,10))
plt.plot(nH,Z,ls='--')
plt.scatter(nH,Z)
plt.hlines(-0.41271532,-7,-3,ls='--',lw=2,color='black')
plt.vlines(-3,-4,-0.41271532,ls='--',lw=2,color='black')
plt.ylabel(r'$\mathbf{log \ [Z/{Z_{\odot}}]}$',labelpad=15,fontsize=30)
# plt.ylabel(ion_label('log \ [Z/{Z_{\odot}}]',''),labelpad=15,fontsize=30)
plt.xlabel(r'$\mathbf{log \ [n_H \ ({cm}^{-3})]}$',labelpad=15,fontsize=30)
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
plt.xlim(-5.25,0.25)
plt.ylim(-1.3,-0.2)
plt.savefig('comp_II_CIE1.png')

plt.show()

quit()

plt.figure()



    
for i,ion in enumerate(ions):
    plt.plot(nH,col_den_ions_dict[ion],ls=ls[i],label=f'{ions_label[i]} ({upper_lim[ion]})')

plt.legend()
plt.ylabel(r'$\mathbf{log \ [N \ {cm}^{-2}]}$',labelpad=15)
plt.xlabel(r'$\mathbf{log \ [n_H \ ({cm}^{-3})]}$',labelpad=15)

plt.show()

quit()


# %%
def physical_param(NHI,NHII,log_T,log_nH):

    log_NH=log10(NHI+NHII)
    l=((NHI+NHII)/(10**log_nH))*(3.2408e-22)
    P_k=10**(log_nH+log_T)

    print(f'NH = {log_NH}    l = {l} kpc     P = {P_k}')


physical_param(1.34896e+14,4.11841e+19,5.29,-3)

# %%
plt.scatter(1,2)
a='O'
b='VI'
c='1036'
plt.text(1,2,f'\\fontsize{{25pt}}{{3em}}\selectfont{{}}$\mathbf{{{a}}}$ \\fontsize{{17pt}}{{3em}}\selectfont{{}}$\mathbf{{{b}}}$ \\fontsize{{25pt}}{{3em}}\selectfont{{}}$\mathbf{{{c}}}$')
# plt.text(1,2,r'{\fontsize{20pt}{3em}\selectfont{}$\mathbf{a}$} {\fontsize{12pt}{3em}\selectfont{}$\mathbf{a}$')
# plt.xlabel(r'{\fontsize{50pt}{3em}\selectfont{}$\mathbf{a}$} {\fontsize{20pt}{3em}\selectfont{}$\mathbf{a}$')
plt.show()


