from operator import indexOf
from astropy.io import fits,ascii
from astropy.table import Table
from numpy import *
from roman import toRoman
from scipy.interpolate import interp2d,interp1d
import pickle
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings("ignore")

ions=['Si+2', 'Si+','C+2', 'C+','O+5']
ions_roman=[f'Si {toRoman(3)}', f'Si {toRoman(2)}',f'C {toRoman(3)}', f'C {toRoman(2)}',f'O {toRoman(6)}']
observations={'Si+2':[12.87,0.08],'Si+':[13.19,0.41], 'C+2':[13.81,0.04],'C+':[14.21,0.39], 'O+5':[13.91,0.04]}

interp_func_dict={}

for i in observations.keys():

    with open(f'Interp_2d_func/{i}_quintic.pkl','rb') as pickle_file:
        f=pickle.load(pickle_file)
    
    interp_func_dict[i]=f

obs_col_den=array([observations[i][0]  for i in ions])
col_den_error=array([observations[i][1]  for i in ions]) 

def chi_sq_func(nH,Z):

    if type(nH)==type(1) or type(nH)==type(1.0) or type(nH)==type(float64(1)):
        nH=[nH]

    else:
        nH=nH

    col_den=array([interp_func_dict[i](nH,Z) for i in observations.keys()])
    chi_sq_exc=zeros(len(nH))
    chi_sq_inc=zeros(len(nH))

    for i in range(len(nH)):

        med_col_den=col_den[:,i]
        chi_sq_exc[i]=sum((((obs_col_den-med_col_den)/col_den_error)**2)[:-1])
        chi_sq_inc[i]=sum((((obs_col_den-med_col_den)/col_den_error)**2))

    return chi_sq_exc,chi_sq_inc


nH_global=arange(-5,0.01,0.01)
Z_global=arange(-3,2.01,0.01)

x=array(list(nH_global)*len(Z_global))
y=zeros(len(x))

chi_sq_exc=zeros((len(Z_global),len(nH_global)))
chi_sq_inc=zeros((len(Z_global),len(nH_global)))

for i,a in enumerate(Z_global):

    chi_sq_exc[i],chi_sq_inc[i]=chi_sq_func(nH_global,a)
    y[len(nH_global)*i:len(nH_global)*(i+1)]=a


def model(chi_sq,name,c):

    chi_sq_arr=chi_sq.flatten()
    ind=argmin(chi_sq_arr)

    nH=round(x[ind],2)
    Z=round(y[ind],2)
    chi_sq_min=round(chi_sq_arr[ind],3)

    if name[:3]=='Exc':
        a=0
    else:
        a=1
    
    chi_sq_err_nH=[chi_sq_func(n,Z)[a] for n in nH_global]
    chi_sq_err_Z=[chi_sq_func(nH,z)[a] for z in Z_global]

    nH_level=[]
    Z_level=[]

    i=indexOf(nH_global,x[ind])

    for j in range(i):
        if chi_sq_err_nH[i-j]-chi_sq_min>1:
            nH_level.append(round(nH_global[i-j],2))
            break

    for j in range(i,len(nH_global)):
        if chi_sq_err_nH[j]-chi_sq_min>1:
            nH_level.append(round(nH_global[j],2))
            break

    i=indexOf(Z_global,y[ind])


    for j in range(i):
        if chi_sq_err_Z[i-j]-chi_sq_min>1:
            Z_level.append(round(Z_global[i-j],2))
            break

    for j in range(i,len(Z_global)):
        if chi_sq_err_Z[j]-chi_sq_min>1:
            Z_level.append(round(Z_global[j],2))
            break

    nH_err=[round(nH-nH_level[0],2),round(nH_level[1]-nH,2)]
    Z_err=[round(Z-Z_level[0],2),round(Z_level[1]-Z,2)]

    print(f'{name} : nH = {nH}  [{nH_err[0]},{nH_err[1]}]   Z = {Z}  [{Z_err[0]},{Z_err[1]}]  chi-sq = {chi_sq_min}')

    col_den=array([interp_func_dict[ion](nH,Z) for ion in observations.keys()])
    xaxis=linspace(1,len(ions_roman),len(ions_roman))

    plt.plot(xaxis,col_den,label=name,ls='--',lw=3,color=c)

    return float(nH),float(Z),nH_err,Z_err


def plot_samples(m,c,n=10,i1=1,i2=1):

    nH_sample=random.normal(m[0],max(m[2])*i1,n)
    Z_sample=random.normal(m[1],max(m[3])*i2,n)

    sample = vstack((nH_sample, Z_sample)).T 
    print(sample)
    xaxis=linspace(1,len(ions_roman),len(ions_roman))

    for s in sample:
        nH,Z=s
        col_den=array([interp_func_dict[i](round(nH,3),round(Z,3)) for i in observations.keys()]) 
        plt.plot(xaxis,col_den,alpha=0.2,color=c)


xaxis=linspace(1,len(ions_roman),len(ions_roman))

plt.figure()

plt.errorbar(xaxis,obs_col_den,c='red',yerr=col_den_error, fmt='o',capsize=3,label='Observed')
m1=model(chi_sq_exc,'Excluding OVI','orange')
m2=model(chi_sq_inc,'Including OVI','green')
plot_samples(m1,'orange')
plot_samples(m2,'green')
plt.xticks(xaxis,ions_roman,fontsize=20)
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







quit()


def error_estimate(nH,Z,case):

    if case=='Exc':
        a=0
    else:
        a=1
    
    chi_sq_min=chi_sq_func(nH,Z)[a]
    chi_sq_err_nH=[chi_sq_func(n,Z)[a] for n in nH_global]
    chi_sq_err_Z=[chi_sq_func(nH,z)[a] for z in Z_global]



    plt.figure()
    plt.title('nH')
    plt.scatter(nH_global,chi_sq_err_nH)
    plt.hlines(chi_sq_min+1,xmin=nH_global[0],xmax=nH_global[-1])

    plt.figure()
    plt.title('Z')
    plt.scatter(Z_global,chi_sq_err_Z)
    plt.hlines(chi_sq_min+1,xmin=Z_global[0],xmax=Z_global[-1])
 











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





































