from astropy.io import fits,ascii
from astropy.table import Table
from numpy import *
from roman import toRoman
from scipy.interpolate import interp2d,interp1d
import pickle
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings("ignore")
plt.style.use('../Python/my_style.mpl')


hdu=fits.open(f'Data/component_III_nH_Z_col_density_param.fits')
data=Table(hdu[1].data)


log_nH=data['log_nH']
log_Z=data['log_Z']

lognH_range=[-5,0]
logZ_range=[-3,2]

param=data['parameters']
ions=['Si+2', 'Si+','C+2', 'C+','O+5']

observations={'Si+2':[12.87,0.08],'Si+':[13.19,0.41], 'C+2':[13.81,0.04],'C+':[14.21,0.39], 'O+5':[13.91,0.04]}


interp_func_dict={}

for i in observations.keys():

    with open(f'Interp_2d_func/{i}_quintic.pkl','rb') as pickle_file:
        f=pickle.load(pickle_file)

    interp_func_dict[i]=f
    

def log_posterior(theta,ions_to_use=['Si+2', 'Si+','C+2', 'C+']):

    lognH,logZ=theta
    
    #prior

    if lognH_range[0]<=lognH<=lognH_range[1] and logZ_range[0]<=logZ<=logZ_range[1]:
         log_prior=0
    
    else:
        log_prior=-inf

    #likelihood
    
    model_col_den=[]
    observed_col_den=[]
    col_den_error=[]

    for i in ions_to_use:

        f=interp_func_dict[i]
        col_den=f(lognH,logZ)
        model_col_den.append(col_den)

        observed_col_den.append(observations[i][0])
        col_den_error.append(observations[i][1])
    
    model_col_den=array(model_col_den)    
    observed_col_den=array(observed_col_den)
    col_den_error=array(col_den_error)

    log_likelihood=-0.5*sum(log(2*pi*(col_den_error**2))+(((observed_col_den-model_col_den)/(col_den_error))**2))

    return log_prior+log_likelihood    # posterior


nH=arange(-5,0,0.01)
Z=arange(-3,2,0.01)

x=zeros(len(nH)*len(Z))
y=zeros(len(nH)*len(Z))
log_post=zeros(len(nH)*len(Z))

k=0
for i,n in enumerate(nH):
    for j,z in enumerate(Z):
        n=round(n,2)
        z=round(z,2)
        x[k]=n
        y[k]=z
        log_post[k]=log_posterior([n,z])
        k+=1


i=argmax(log_post)

print(f'nH = {x[i]}  Z = {y[i]}')

fig=plt.figure()
ax=plt.axes(projection ='3d')

ax.scatter(x,y,log_post)
ax.set_xlabel('nH')
ax.set_ylabel('Z')
ax.set_zlabel('log_post')

plt.show()






quit()

ions=['Si+2', 'Si+','C+2', 'C+','O+5']
observations={'Si+2':[12.87,0.08],'Si+':[13.19,0.41], 'C+2':[13.81,0.04],'C+':[14.21,0.39], 'O+5':[13.91,0.04]}

i='O+5'

with open(f'Interp_2d_func/{i}_quintic.pkl','rb') as pickle_file:
    f=pickle.load(pickle_file)



obs_col_den=array([observations[i][0]  for i in ions])
col_den_error=array([observations[i][1]  for i in ions]) 


nH=arange(-5,0,0.01)
Z=arange(-3,2,0.01)


x=array(list(nH)*len(Z))
y=zeros(len(x))

col_den=zeros((len(Z),len(nH)))

for i,z in enumerate(Z):

    col_den[i]=f(nH,z)
    y[len(nH)*i:len(nH)*(i+1)]=z


z=col_den.flatten()

fig=plt.figure()
ax=plt.axes(projection ='3d')

ax.scatter(x,y,z)
ax.set_xlabel('nH')
ax.set_ylabel('Z')
ax.set_zlabel('col_den')

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





































