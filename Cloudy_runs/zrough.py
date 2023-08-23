from astropy.io import fits,ascii
from astropy.table import Table
from numpy import *
from scipy.interpolate import interp2d,interp1d
import pickle
import matplotlib.pyplot as plt




'plotting saved functions'

def plot_saved_func():
 
    i='O+5'


    with open(f'Interp_2d_func/{i}_quintic.pkl','rb') as pickle_file:
            f=pickle.load(pickle_file)

    z=linspace(-3,2,100)
    nH=arange(-5,0.1,1)
    marker=['D','+','x','s','o','^','*']

    for i,n in enumerate(nH):
        plt.scatter(z,f(n,z),label=f'{n}',marker=marker[i])

    plt.legend()
    plt.show()
    quit()


# plt.style.use('../Python/my_style.mpl')

hdu=fits.open(f'Data/component_III_nH_Z_col_density_param.fits')
data_2d=Table(hdu[1].data)

hdu=fits.open(f'Data/component_III_nH_col_density_param.fits')
data_1d=Table(hdu[1].data)




log_nH_1d=data_1d['log_nH']

log_nH_2d=data_2d['log_nH']
log_Z_2d=data_2d['log_Z']


lognH_range=[-5,0]
logZ_range=[-3,0]
log_Z_ref=-1

param=data_2d['parameters']
ions=['Si+', 'Si+2','C+', 'C+2','O+5']

observations={'Si+':[13.19,0.41], 'Si+2':[12.87,0.08],'C+':[14.21,0.39], 'C+2':[13.81,0.04],'O+5':[13.91,0.04]}


for i in observations.keys():

    plt.figure()
    plt.scatter(log_nH_1d,log10(data_1d[i]))
    plt.title(i)

# plt.show() hdhldh
quit()

def save_func(func,name):

    with open(f'Interp_2d_func/{name}_quintic.pkl','wb') as pickle_file:
        pickle.dump(func,pickle_file)

def interp_func(ion,kind='cubic'):

    ion_col_den_1d=data_1d[ion].value

    for i,j in enumerate(ion_col_den_1d):
        if j==0:
            break

    log_nH_1df=log_nH_1d[:i]
    log_col_den=log10(ion_col_den_1d[:i])
    f_1d=interp1d(log_nH_1df,log_col_den,fill_value='extrapolate',kind='quadratic')

    ion_col_den_2d=data_2d[ion].value
    log_col_den=zeros(len(data_2d))

    for i,j in enumerate(ion_col_den_2d):

        if j!=0:
            log_col_den[i]=log10(j)

        else:
            log_col_den[i]=f_1d(log_nH_2d[i])+log_Z_2d[i]-log_Z_ref
            

    f_2d=interp2d(log_nH_2d,log_Z_2d,log_col_den,kind=kind)
    save_func(f_2d,ion)

    return f_2d

def interp_func_1d():

    func_dict={}
    # log_nH=data['log_nH']

    for ion in observations.keys():

        ion_col_den=data_1d[ion].value

        for i,j in enumerate(ion_col_den):
            if j==0:
                break

        log_nH_f=log_nH_1d[:i]
        log_col_den=log10(ion_col_den[:i])
        f=interp1d(log_nH_f,log_col_den,fill_value='extrapolate',kind='quadratic')
        func_dict[ion]=f

    return func_dict

interp_func_1d_dict=interp_func_1d()

def plot_1d_interp(ion,nH,Z):

    log_col_den=interp_func_1d_dict[ion](nH)+Z+1

    plt.scatter(Z,log_col_den,label=f'{nH} 1d')






def col_vs_z(ion,nH):

    mask=log_nH_2d==float(nH)
    data1=data_2d[mask]
    Z=data1['log_Z']
    log_col_den=log10(data1[ion])

    plt.scatter(Z,log_col_den,label=f'{nH}',marker='D')  
    
i='O+5'

# with open(f'Interp_2d_func/{i}_quintic.pkl','rb') as pickle_file:
#     f=pickle.load(pickle_file)

z=linspace(-3,2,100)
nH=arange(-5,0,1)

plt.figure()

for n in nH:
    # plt.scatter(z,f(n,z),label=f'{n}')
    col_vs_z(i,n)
    plot_1d_interp(i,n,z)


plt.legend()
# plt.title(i)

plt.show()





'interp1 function'

# hdu=fits.open(f'Data/component_III_nH_col_density_param.fits')
# data=Table(hdu[1].data)

# nH=data['log_nH']

# oVI=data['C+']

# def interp_func(ion_col_den):

#     for i,j in enumerate(ion_col_den):
#         if j==0:
#             break

#     log_nH=nH[:i]
#     log_col_den=log10(ion_col_den[:i])

#     f=interp1d(log_nH,log_col_den,fill_value='extrapolate',kind='quadratic')

#     return f
























quit()

# data=loadtxt('text.txt',dtype=str)

# a=data[:,1]

# v=[float(x[:-2]) for x in a]

# print(v)


'observed and predicted'


# plt.style.use('../Python/my_style.mpl')


# ions=[f'Si {toRoman(2)}', f'Si {toRoman(3)}',f'C {toRoman(2)}', f'C {toRoman(3)}', f'O {toRoman(6)}']
# x=linspace(1,len(ions),len(ions))


# obs_col_density=array([12.53,12.69,13.49,13.65,13.84])
# col_density_error=array([0.07,0.06,0.05,0.03,0.03])

# predicted_all_ions=log10([4.43676e+11, 6.40882e+11, 1.18221e+14, 1.23953e+14, 4.68142e+13])
# predicted_without_OVI=log10([3.29517e+12, 2.07164e+12,1.02575e+14, 3.75366e+13, 4.88502e+11])

# fig=plt.figure(figsize=(11,11))

# plt.errorbar(x,obs_col_density,c='red',yerr=col_density_error, fmt='o',capsize=3,label='Observed')
# plt.plot(x,predicted_all_ions,label=r'All ions $(T={10}^{4.49}\ K)$',ls='--')
# plt.plot(x,predicted_without_OVI,label=r'Excluding OVI $(T={10}^{4.38}\ K)$',ls='--')
# plt.xticks(x,ions,fontsize=20)
# plt.yticks(fontsize=20)
# plt.ylabel(r'$\mathbf{log \ (N \ {cm}^{-2})}$',labelpad=15)
# plt.xlabel(r'$\mathbf{Ions}$',labelpad=15)
# plt.legend(loc='upper left')
# plt.savefig('Files_n_figures/Observed_and_predicted.png')
# plt.show()

# quit()

