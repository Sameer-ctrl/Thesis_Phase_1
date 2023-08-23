from cProfile import label
from traceback import print_tb
from astropy.io import fits,ascii
from astropy.table import Table
from numpy import *
from roman import toRoman
from scipy.interpolate import interp2d,interp1d
import pickle
import matplotlib.pyplot as plt

ions=['Si+2', 'Si+','C+2', 'C+','O+5']
observations={'Si+2':[12.87,0.08],'Si+':[13.19,0.41], 'C+2':[13.81,0.04],'C+':[14.21,0.39], 'O+5':[13.91,0.04]}


interp_func_dict={}

for i in observations.keys():

    with open(f'Interp_2d_func/{i}_quintic.pkl','rb') as pickle_file:
        f=pickle.load(pickle_file)

    interp_func_dict[i]=f


obs_col_den=array([observations[i][0]  for i in ions])
col_den_error=array([observations[i][1]  for i in ions]) 

def chi_sq(nH,Z):

    if type(nH)==type(1) or type(nH)==type(1.0):
        nH=[nH]

    else:
        nH=nH

    median_col_den_exc_OVI=array([interp_func_dict[i](nH,Z) for i in observations.keys()])
    chi_sq_exc=zeros(len(nH))

    for i in range(len(nH)):

        med_col_den=median_col_den_exc_OVI[:,i]
        chi_sq_exc[i]=sum((((obs_col_den-med_col_den)/col_den_error)**2)[:-1])

    return chi_sq_exc


nH=arange(-5,0,0.01)
Z=arange(-3,2,0.01)

# nH=linspace(-5,0,100)
# Z=linspace(-3,2,100)

x=array(list(nH)*len(Z))
y=zeros(len(x))

chi_sq_arr=zeros((len(Z),len(nH)))

for i,z in enumerate(Z):

    chi_sq_arr[i]=chi_sq(nH,z)
    y[len(nH)*i:len(nH)*(i+1)]=z


z=chi_sq_arr.flatten()
i=argmin(z)

nH=round(x[i],2)
Z=round(y[i],2)
print('\n ----------------------------------- \n')
print(f'nH = {nH}  Z = {Z}  chi-sq = {round(z[i],3)}')
print('\n ----------------------------------- \n')

ions_all=[f'Si {toRoman(3)}', f'Si {toRoman(2)}',f'C {toRoman(3)}', f'C {toRoman(2)}',f'O {toRoman(6)}']
median_col_den_exc_OVI=array([interp_func_dict[i](nH,Z) for i in observations.keys()])

obs_col_den=array([observations[i][0]  for i in ions])
col_den_error=array([observations[i][1]  for i in ions])


x=linspace(1,len(ions_all),len(ions_all))

plt.errorbar(x,obs_col_den,c='red',yerr=col_den_error, fmt='o',capsize=3,label='Observed')
plt.plot(x,median_col_den_exc_OVI,label='median solution (excluding O VI)',ls='--',lw=3,color='orange')
plt.xticks(x,ions_all,fontsize=20)
plt.yticks(fontsize=20)
plt.ylabel(r'$\mathbf{log \ (N \ {cm}^{-2})}$',labelpad=15)
plt.xlabel(r'$\mathbf{Ions}$',labelpad=15)
plt.legend()
# plt.savefig('Files_n_figures/Observed_and_predicted.png')
plt.show()






# fig=plt.figure()
# ax=plt.axes(projection ='3d')

# ax.scatter(x,y,z)
# ax.set_xlabel('nH')
# ax.set_ylabel('Z')
# ax.set_zlabel('chi_square')



plt.show()


















quit()

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

# ind=[]

# for i,row in enumerate(data_2d):

#     if row['H']==0 or row['H']<10**16:
#        ind.append(i) 

# data_2d.remove_rows([ind])

log_nH_1d=data_1d['log_nH']

log_nH_2d=data_2d['log_nH']
log_Z_2d=data_2d['log_Z']


lognH_range=[-5,0]
logZ_range=[-3,0]
log_Z_ref=-1

param=data_2d['parameters']
ions=['Si+', 'Si+2','C+', 'C+2','O+5']

observations={'Si+':[13.19,0.41], 'Si+2':[12.87,0.08],'C+':[14.21,0.39], 'C+2':[13.81,0.04],'O+5':[13.91,0.04]}

def save_func(func,name):

    with open(f'Data/{name}.pkl','wb') as pickle_file:
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


def col_vs_z(ion,nH):

    mask=log_nH_2d==float(nH)
    data1=data_2d[mask]
    Z=data1['log_Z']
    log_col_den=log10(data1[ion])

    plt.scatter(Z,log_col_den,label=f'{nH}',marker='D')  
    
i='O+5'

with open(f'Data/{i}.pkl','rb') as pickle_file:
    f=pickle.load(pickle_file)

z=linspace(-3,2,100)
nH=arange(-5,0.1,1)

plt.figure()

for n in nH:
    plt.scatter(z,f(n,z),label=f'{n}')
    col_vs_z(i,n)

plt.legend()
plt.title(i)

plt.show()





































