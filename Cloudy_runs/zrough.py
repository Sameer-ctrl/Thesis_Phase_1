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


'------------------------------------'
'Checking interpolation function'


hdu=fits.open('Data/component_II_T_col_density_param.fits')
data=Table(hdu[1].data)

T=data['log_T']
col_den_OVI=log10(data['O+5'])

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









'------------------------------------'
'Likelihood'

hdu=fits.open(f'Data/component_III_nH_Z0_col_density_param.fits')
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

'------------------------------------'

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











# hdu=fits.open(f'Data/component_III_nH_Z_col_density_param.fits')
# data=Table(hdu[1].data)

# ind=[]

# for i,row in enumerate(data):

#     if row['H']==0 or row['H']<10**16 or row['log_nH']>0 or row['log_Z']>0:
#        ind.append(i) 

# data.remove_rows([ind])



# H=data['H']

# plt.hist(log10(H),bins='auto')



# plt.scatter(nH,Z)
# plt.show()








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


























