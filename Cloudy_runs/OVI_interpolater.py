from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
from numpy import *
from scipy.interpolate import interp2d
import pickle


hdu=fits.open(f'pg0003/z=0.347579/nH_Z/component_III_nH_Z_col_density_param.fits')
data_2d=Table(hdu[1].data)

hdu=fits.open(f'pg0003/z=0.347579/nH_Z=-1/component_III_nH_col_density_param.fits')
data_1d=Table(hdu[1].data)

log_nH_2d=data_2d['log_nH']
log_Z_2d=data_2d['log_Z']

ions=['Si+', 'Si+2','C+', 'C+2','O+5']

for i in ions:

    col_den_OVI=log10(data_2d[i])
    plt.figure(figsize=(10,10))
    plt.scatter(log_nH_2d,log_Z_2d,c=col_den_OVI)
    plt.colorbar()
    plt.xlim(-5,0)
    plt.ylim(-3,2)
    plt.title(i)
    plt.savefig(f'{i}_obs.png')


# plt.show()



quit()



hdu=fits.open('Data/col_den_OVI.fits')
data=Table(hdu[1].data)

log_nH=data['log_nH']
log_Z=data['log_Z']
log_col_den=data['O+5']

kind='quintic'

# f=interp2d(log_nH,log_Z,log_col_den,kind=kind)


# with open(f'O+5_{kind}.pkl','wb') as pickle_file:
#         pickle.dump(f,pickle_file)

def col_vs_z(nH):

    mask=log_nH==float(nH)
    data1=data[mask]
    Z=data1['log_Z']
    log_col_den=log10(data1['O+5'])

    plt.scatter(Z,log_col_den,label=f'{nH}',marker='D')  


i='O+5'

with open(f'{i}_{kind}.pkl','rb') as pickle_file:
    f=pickle.load(pickle_file)

with open(f'{i}_{kind}1.pkl','rb') as pickle_file:
    f1=pickle.load(pickle_file)

z=linspace(-3,2,100)
nH=arange(-5,0.1,0.5)

# plt.figure(figsize=(16,9))

# for n in nH:
#     n=round(n,2)
#     plt.scatter(z,f(n,z),label=f'{n}')
#     # col_vs_z(n)

# plt.legend()
# plt.title(i)
# plt.show() 

nH=linspace(-5,0,500)
Z=linspace(-3,2,500)

x=array(list(nH)*len(Z))
y=zeros(len(x))
col_den=zeros((len(Z),len(nH)))
col_den1=zeros((len(Z),len(nH)))


for i,z in enumerate(Z):

    col_den[i]=f(nH,z)
    col_den1[i]=f1(nH,z)
    y[len(nH)*i:len(nH)*(i+1)]=z


col_den=col_den.flatten()
col_den1=col_den1.flatten()


fig=plt.figure()
ax=plt.axes(projection ='3d')   

ax.scatter(log_nH,log_Z,log_col_den,label='interp1')
ax.scatter(x,y,col_den,label='interpolated')
ax.scatter(x,y,col_den1,label='interpolated2')
ax.set_xlabel('nH')
ax.set_ylabel('Z')
ax.set_zlabel('col den')

plt.legend()

plt.show()