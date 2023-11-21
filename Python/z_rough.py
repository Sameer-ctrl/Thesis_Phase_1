import os
from astropy.io import fits,ascii
from astropy.table import Table,Column
import matplotlib.pyplot as plt
from numpy import *
from scipy.integrate import simpson

files=os.listdir('Files_n_figures/sys_plots')
files=sorted(files)

for file in files:
    file_split=file.split('_')

    print(f'\\begin{{figure}} \n  \centering  \n  \hspace*{{-21mm}}\n    \captionsetup{{oneside,margin={{0cm,21mm}}}}\n    \includegraphics[width=\linewidth]{{Figures//system-plots/{file}}} \n  \caption{{System plot of the BLA candidate towards LOS of {file_split[0]} at $z_{{abs}}=${file_split[1][2:]}}} \n\end{{figure}}')
    print('\n\n')

# \begin{figure}
#     \centering
#     \includegraphics[width=1\linewidth]{Figures//system-plots/1ES1553+113_z=0.187764_sys_plot.png}
#     \caption{Enter Caption}
# \end{figure}



quit()

plt.style.use('my_style.mpl')

def ion_label(ion):

    for i,s in enumerate(ion):
        if s.isupper():
            if i==0:
                pass
            else:
                break

    atom=ion[:i]
    n=ion[i:]

    return f'{{\\fontsize{{25pt}}{{3em}}\selectfont{{}}$\mathbf{{{atom}}}$}} {{\\fontsize{{17pt}}{{3em}}\selectfont{{}}$\mathbf{{{n}}}$}}'


total_ions=['CIII', 'NV', 'OVI', 'CIV', 'SiIV', 'SiIII', 'OVI', 'CIV', 'SiIII', 'NV', 'CIV', 'SiIV', 'NV', 'OVI', 'SiIII', 'CIII', 'OVI', 'SiIII', 'NV', 'CIV', 'SiIII', 'NV', 'CIV', 'SiIV', 'SiIII', 'CIII', 'NV', 'OVI', 'OIII', 'CIII', 'OVI', 'OIII', 'CII', 'NV', 'FeII', 'SiIII', 'AlII', 'CIV', 'SiIV', 'OI', 'SiII', 'CII*', 'NII', 'NV', 'CII', 'SiIII', 'OVI', 'CIV', 'SiIV', 'PII', 'SiII', 'CIII', 'CII', 'SiIII', 'OVI', 'SiIV', 'SiII', 'CII', 'SiII', 'OI', 'CIII', 'OVI', 'SiIII', 'CIV', 'SiIV', 'SiIII', 'CIII', 'OVI', 'OIII', 'SiIII', 'CIV', 'SiIV', 'SiIII', 'SiIV', 'OVI', 'CIV', 'SiIII', 'NII', 'CII', 'FeII', 'CIV', 'SiIV', 'OI', 'SiII', 'NII', 'NV', 'NIII', 'CIII', 'CII', 'SiIII', 'OVI', 'SiIV', 'OI', 'SiII', 'NV', 'OVI', 'SiIII', 'CIII', 'OVI', 'SiIII', 'CII', 'SiIII', 'SiII', 'CIV', 'SiIV', 'SiIII', 'CII', 'SiIII', 'SiIV', 'OVI', 'CIV', 'CIII', 'CII', 'NIII', 'SiIII', 'OVI', 'OI', 'SiII', 'CII', 'FeII', 'SiII']


hist_plot=plt.hist(total_ions,bins=len(set(total_ions)),rwidth=0.5,color='#a1c9f4')
x=array([ 0 ,  0.9375,  1.875 ,  2.8125,  3.75  ,  4.6875,  5.625 , 6.5625,  7.5   ,  8.4375,  9.375 , 10.3125, 11.25  , 12.1875, 13.125 , 14.0625, 15])# ions=plt.xticks()
ions_ticks= [  'CIII',  'NV',  'OVI',  'CIV',  'SiIV',  'SiIII',  'OIII',  'CII', 'FeII', 'AlII','OI',  'SiII', 'CII*', 'NII',  'PII', 'NIII']
ions_ticks=[ion_label(x) for x in ions_ticks]
plt.xticks(x[:-1]+(0.9375/2),ions_ticks)
plt.yticks([0,4,8,12,16,20,24],fontsize=25)
plt.xlabel(r'$\mathbf{Metal \ Ions}$',labelpad=20,fontsize=30)
plt.ylabel(r'$\mathbf{n}$',labelpad=20,fontsize=30)
plt.show()

# a=f'{{\\fontsize{{25pt}}{{3em}}\selectfont{{}}$\mathbf{{{line_name[0]}}}$}} {{\\fontsize{{17pt}}{{3em}}\selectfont{{}}$\mathbf{{{line_name[1]}}}$}}'





quit()

plt.style.use('my_style.mpl')


b=linspace(20,100,5,dtype=int)
# N=[13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18,18.5,19]
N=[13,13.5,14,15,16.5,17.5,18.5]
a=1.25
plt.figure(figsize=(15*a,9*a),dpi=300)

# for val in N:
for val in b:
    # file=f'voigt_mock/N_{val}.txt'
    file=f'voigt_mock/b_{val}_1.txt'

    data=loadtxt(file,comments='!')

    cen_wave_rest=1215.6701
    cen_wave_obs=cen_wave_rest*(1)

    wave=data[:,0]
    cont=data[:,3]
    v=3*(10**5)*((wave**2-(cen_wave_obs**2))/(wave**2+(cen_wave_obs**2)))

    plt.plot(v,cont,label=f'{{$\mathbf{{b=}} \ \mathbf{{{val}}} \ \mathbf{{km \ s^{{-1}}}}$}}',lw=2)
    # plt.plot(v*(3/5),cont,label=f'{{$\mathbf{{N=}} \ \mathbf{{{val}}} \ \mathbf{{cm^{{-2}}}} $}}',lw=2)

    # v2=3*(10**5)*((wave_spec**2-(cen_wave_obs**2))/(wave_spec**2+(cen_wave_obs**2)))
# \mathbf{10^{13} \ cm^{-2}}
plt.legend(loc='lower right',fontsize=22*a)
plt.xlim(-300,300)
# plt.ylim(bottom=0.65)
plt.xticks(fontsize=20*a)
plt.yticks(fontsize=20*a)
plt.xlabel(r'$\mathbf{V} \ \mathbf{(km \ \ s^{-1})}$',labelpad=20,fontsize=30*a)
plt.ylabel(r'$\mathbf{Normalized \ Flux} $',labelpad=20,fontsize=30*a)
# plt.text(-270,0.185,r'$\mathbf{{b=}} \ \mathbf{{{50}}} \ \mathbf{{km \ s^{{-1}}}}$',fontsize=25*a)
plt.text(-250,0.725,r'$\mathbf{{N=}} \ \mathbf{{10^{13}}} \ \mathbf{{cm^{{-2}}}}$',fontsize=25*a)
plt.savefig('Voigt-b.png')
# plt.show()


# wave=linspace(1200,1230,10000)
# flux=ones(10000)
# err=ones(10000)*0.1

# wave_col=Column(name='WAVE',data=wave)
# flux_col=Column(name='FLUX',data=flux)
# err_col=Column(name='ERROR',data=err)

# tab=Table()

# tab.add_columns([wave_col,flux_col,err_col])
# tab.write(f'Mock.asc', format='ascii', overwrite=True)





quit()

data=loadtxt('Data/rest_wave.txt',dtype=str)

ion=data[:,1][13:26]
wave=data[:,0].astype(float)[13:26]

for i in range(len(ion)):
    print(f'{ion[i]}  {1347.4/wave[i]-(1):f}')


quit()


file_path='Data/VPfit_fits_rebinned/Metals_HI'

data=loadtxt(f'{file_path}/HI_1215.txt',comments='!')
wave=data[:,0]
cont=data[:,3]

data1=loadtxt(f'{file_path}/HI_1215_1.txt',comments='!')
wave1=data1[:,0]
cont1=data1[:,3]

data2=loadtxt(f'{file_path}/HI_1215_2.txt',comments='!')
wave2=data2[:,0]
cont2=data2[:,3]

data3=loadtxt(f'{file_path}/HI_1215_3.txt',comments='!')
wave3=data3[:,0]
cont3=data3[:,3]


plt.plot(wave,cont)
plt.plot(wave1,cont1,ls='--')
plt.plot(wave2,cont2,ls='--')
plt.plot(wave3,cont3,ls='--')
plt.plot(wave,(cont1*cont2*cont3),ls='-.')

plt.show()






quit()

z_abs=0.347

z=array([0.347,0.352484])

# v=3e5*(((1+z)**2-1)/((1+z)**2+1))

# print(v[0]-v[1])
z1=0.347
z2=0.34835

del_v=3e5*((1+z1)**2-(1+z2)**2)/((1+z1)**2+(1+z2)**2)
print(del_v)


quit()

# 'comparison of different continuum, binning'

# low=loadtxt('Data/fit_param_binned_low.txt',dtype=str)
# high=loadtxt('Data/fit_param_binned_high.txt',dtype=str)
# mean_cont=loadtxt('Data/fit_param_binned_mean.txt',dtype=str)
# unbinned=loadtxt('Data/fit_param_unbinned.txt',dtype=str)

# class fits_param():

#     def __init__(self,data):
        
#         self.lines=data[1:,0]
#         self.z=data[1:,1].astype(float)
#         self.b=data[1:,3].astype(float)
#         self.logN=data[1:,5].astype(float)
#         self.err_z=data[1:,2].astype(float)
#         self.err_b=data[1:,4].astype(float)
#         self.err_logN=data[1:,6].astype(float)
#         self.chi_sq=data[1:,7].astype(float)


# fits_param_low=fits_param(low)
# fits_param_high=fits_param(high)
# fits_param_mean=fits_param(mean_cont)

# err_ml=fits_param_mean.logN-fits_param_low.logN
# err_mh=fits_param_mean.logN-fits_param_high.logN

# err_ml=[abs(round(x,2)) for x in err_ml]
# err_mh=[abs(round(x,2)) for x in err_mh]

# print(err_ml)
# print(err_mh)



#     def plot_param(self,param,spec):

#         x=linspace(1,len(self.lines),len(self.lines))

#         if param=='z':
#             y=self.z
#             err=self.err_z
        
#         elif param=='b':
#             y=self.b
#             err=self.err_b
        
#         elif param=='logN':
#             y=self.logN
#             err=self.err_logN
        
#         elif param=='chi_sq':
#             y=self.chi_sq
            

#         plt.scatter(x,y,label=spec,s=50)
#         plt.plot(x,y,ls='--')
#         plt.ylabel(f'{param}',labelpad=15)
#         plt.xlabel('line',labelpad=15)
#         plt.xticks(x,self.lines)


# fits_param_low=fits_param(low)
# fits_param_high=fits_param(high)
# fits_param_mean=fits_param(mean_cont)




# fits_param_unbinned=fits_param(unbinned)

# plt.figure()

# fits_param_mean.plot_param('chi_sq','rebinned')
# # fits_param_low.plot_param('z','lower continuum')
# # fits_param_high.plot_param('z','upper continuum')
# fits_param_unbinned.plot_param('chi_sq','oversampled')
# plt.title(r'${\chi}^{2}$ values for rebinned and oversampled spectrum')

# plt.legend()
# plt.show()


# quit()


'bin v/s unbinned spectrum plot'


# file='PG0003+158_rebinned.fits'

# data_a=loadtxt('Data/spec_PG0003+158_v3.dat')

# wave_a=data_a['WAVE']
# flux_a=data_a['FLUX']

# hdu_org=fits.open('Data/PG0003+158.fits')
# data_org=Table(hdu_org[1].data)

# wave_org=data_org['WAVE'][data_org['WAVE']>=1132.7]
# flux_org=data_org['FLUX'][data_org['WAVE']>=1132.7]


# hdu=fits.open(f'Data/{file}')
# data=Table(hdu[1].data)

# wave=data['WAVE']
# flux=data['FLUX']
# cont=data['CONT_FLUX']
# err=data['ERROR']


# plt.step(wave_org,flux_org,label='Unbinned spectrum')
# plt.step(wave,flux,label='My rebinned spectrum')
# plt.step(wave_a,flux_a,label='Your rebinned spectrum')
# # plt.plot(wave,cont,label='mean')
# # plt.plot(wave,cont*0.97,label='lower',ls='--')
# # plt.plot(wave,cont*1.03,label='upper',ls='--')
# plt.ylabel('Flux')
# plt.xlabel('Wavelength')
# plt.legend()
# plt.show()


'equivalent width'


# hdu_unbin=fits.open('Data/PG0003+158_unbinned.fits')
# hdu_bin=fits.open('Data/PG0003+158_rebinned.fits')

# class Spectrum():

#     def __init__(self,hdu):

#         data=Table(hdu[1].data)

#         self.wave=data['WAVE']
#         self.flux=data['FLUX']
#         self.cont=data['CONT_FLUX']
#         self.err=data['ERROR']
    
#     def eqw(self,ions):
        
#         cont=self.cont
#         wave=self.wave
#         flux=self.flux

#         eq_w=[]

#         for i in ions:

#             wave_slice=[]
#             flux_slice=[]
#             cont_slice=[]

#             for j in range(len(wave)):
#                 if ions[i][0] <= wave[j] <= ions[i][1]:
#                     wave_slice.append(wave[j])
#                     flux_slice.append(flux[j])
#                     cont_slice.append(cont[j])

#             integral=simpson(flux_slice,wave_slice)
#             w=(ions[i][1]-ions[i][0])-(integral/mean(cont_slice))

#             eq_w.append(w)

#         return eq_w



# ions={'Ha':[1636.55,1639.225],'Hb':[1381.1,1383.264],'OVI_1031':[1389.856,1391.254],'OVI_1037':[1397.816,1398.877],'CII':[1396.727,1396.995],'CIII':[1310.485,1311.182],'SiII':[1698.778,1699.039],'SiIII':[1626.068,1626.403]}


# spec_unbin=Spectrum(hdu_unbin)
# spec_bin=Spectrum(hdu_bin)


# plt.step(spec_bin.wave,spec_bin.flux)
# plt.show()


# quit()
# eq_w_unbin=array(spec_unbin.eqw(ions))
# eq_w_bin=array(spec_bin.eqw(ions))

# diff=eq_w_bin-eq_w_unbin

# x=linspace(1,len(ions),len(ions))

# plt.plot(x,eq_w_unbin,ls='--')
# plt.plot(x,eq_w_bin,ls='--')
# plt.scatter(x,eq_w_unbin,label='unbinned')
# plt.scatter(x,eq_w_bin,label='rebinned')
# plt.xticks(x,ions.keys())
# plt.xlabel('Lines',labelpad=15)
# plt.ylabel('Equivalent Width',labelpad=15)
# plt.legend()
# plt.show()

# plt.step(spec_unbin.wave,spec_unbin.flux)
# plt.plot(spec_unbin.wave,spec_unbin.cont)
# plt.show()


'Velocity separation of components'

# data=loadtxt('wavelengths_fit.txt',dtype=str)
# data_rest=loadtxt('Data/rest_wave.txt',dtype=str)

# ion=data_rest[:,1]
# wave_rest=data_rest[:,0].astype(float)

# rest_wave={}

# for i in range(len(ion)):
#     rest_wave.update({ion[i]:wave_rest[i]})

# line=data[1:,0]
# wave=data[1:,1].astype(float)
# v=[]

# z_abs=(wave[3]-rest_wave['OVI_1031'])/rest_wave['OVI_1031']

# for i in range(len(line)):

#     cen_wave_obs=(1+z_abs)*rest_wave[line[i]]
#     a=((wave[i]**2-(cen_wave_obs**2))/(wave[i]**2+(cen_wave_obs**2)))
#     v.append(3e5*(a))

# print(line)
# print(v)


'temperature using 0 VI and HI'

# bHI=62.49162
# bOVI=29.63435

# T=(16/(15*(0.129**2)))*(bHI**2-(bOVI**2))
# print(log10(T))


'difference in paramaters due to continuum'




