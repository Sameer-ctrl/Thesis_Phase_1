from astropy.io import fits
from astropy.table import Table,Column
import matplotlib.pyplot as plt
from numpy import *
from scipy.interpolate import interp2d,interp1d
import pickle

def extrapolate_OVI():

    hdu=fits.open(f'pg0003/z=0.347579/nH_Z/component_III_nH_Z_col_density_param.fits')
    data_2d=Table(hdu[1].data)

    hdu=fits.open(f'pg0003/z=0.347579/nH_Z=-1/component_III_nH_col_density_param.fits')
    data_1d=Table(hdu[1].data)

    mask=data_2d['O+5']>0

    log_nH_2d=data_2d['log_nH'][mask]
    nH_2d_unique=unique(data_2d['log_nH'])

    log_Z_2d=data_2d['log_Z'][mask]
    col_den_OVI=log10(data_2d['O+5'][mask])


    table=Table()
    logZ_table=[]
    lognH_table=[]
    col_den_table=[]

    Z=unique(log_Z_2d)

    for z in Z:

        mask=log_Z_2d==z

        if z>=0.64 and z<=0.88 :
            nH=log_nH_2d[mask][:-1]
            col_den_OVI_masked=col_den_OVI[mask][:-1]
        
        else:
            nH=log_nH_2d[mask]
            col_den_OVI_masked=col_den_OVI[mask]
        

        f=interp1d(nH,col_den_OVI_masked,fill_value='extrapolate',kind='quadratic')

        logZ=z*ones(len(nH_2d_unique))
        col_den=f(nH_2d_unique)
        
        logZ_table.append(logZ)
        lognH_table.append(nH_2d_unique)
        col_den_table.append(col_den)



    Z_col=Column(name='log_Z',data=array(logZ_table).flatten())
    nH_col=Column(name='log_nH',data=array(lognH_table).flatten())
    col_den_col=Column(name='O+5',data=array(col_den_table).flatten())

    table.add_columns([nH_col,Z_col,col_den_col])
    table.write('pg0003/z=0.347579/nH_Z/col_den_OVI_extrapolated.fits')



    
def interpolate_OVI(kind='quintic'):

    hdu=fits.open('pg0003/z=0.347579/nH_Z/col_den_OVI_extrapolated.fits')
    data=Table(hdu[1].data)

    log_nH=data['log_nH']
    log_Z=data['log_Z']
    log_col_den=data['O+5']


    f=interp2d(log_nH,log_Z,log_col_den,kind=kind)

    with open(f'Interp_2d_func_new/O+5_{kind}.pkl','wb') as pickle_file:
            pickle.dump(f,pickle_file)


# extrapolate_OVI()
interpolate_OVI(kind='cubic')

quit()


# plt.scatter(nH_col,Z_col,c=col_den_col)
# plt.colorbar()
# plt.show()

# fig=plt.figure()
# ax=plt.axes(projection ='3d')   

# ax.scatter(nH_col,Z_col,col_den_col,label='interp1')
# ax.scatter(x,y,col_den,label='interpolated')
# ax.scatter(x,y,col_den1,label='interpolated2')
# ax.set_xlabel('nH')
# ax.set_ylabel('Z')
# ax.set_zlabel('col den')
# plt.legend()

# plt.show()




hdu=fits.open('pg0003/z=0.347579/nH_Z/col_den_OVI_extrapolated.fits')
data=Table(hdu[1].data)

log_nH=data['log_nH']
log_Z=data['log_Z']
log_col_den=data['O+5']

kind='quintic'


def col_vs_z(nH):

    mask=log_nH==float(nH)
    data1=data[mask]
    Z=data1['log_Z']
    log_col_den=log10(data1['O+5'])

    plt.scatter(Z,log_col_den,label=f'{nH}',marker='D')  


i='O+5'

with open(f'Interp_2d_func_new/O+5_quintic.pkl','rb') as pickle_file:
    f=pickle.load(pickle_file)


# z=linspace(-3,2,100)
# nH=arange(-5,0.1,0.5)

# plt.figure(figsize=(16,9))

# for n in nH:
#     n=round(n,2)
#     plt.scatter(z,f(n,z),label=f'{n}')
#     # col_vs_z(n)

# plt.legend()
# plt.title(i)
# plt.show() 

nH=linspace(-5,0,100)
Z=linspace(-3,2,100)

x=array(list(nH)*len(Z))
y=zeros(len(x))
col_den=zeros((len(Z),len(nH)))
col_den1=zeros((len(Z),len(nH)))


for i,z in enumerate(Z):

    col_den[i]=f(nH,z)
    # col_den1[i]=f1(nH,z)
    y[len(nH)*i:len(nH)*(i+1)]=z


col_den=col_den.flatten()
col_den1=col_den1.flatten()


fig=plt.figure()
ax=plt.axes(projection ='3d')   

ax.scatter(log_nH,log_Z,log_col_den,label='interp1')
# ax.scatter(x,y,col_den,label='interpolated')
# ax.scatter(x,y,col_den1,label='interpolated2')
ax.set_xlabel('nH')
ax.set_ylabel('Z')
ax.set_zlabel('col den')

plt.legend()

plt.show()