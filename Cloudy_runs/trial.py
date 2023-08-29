from operator import indexOf
from astropy.io import fits
import matplotlib.pyplot as plt
from numpy import *
from astropy.table import Table
from scipy.interpolate import interp1d


hdu=fits.open('Data/component_III_nH_Z_col_density_param.fits')
data=Table(hdu[1].data)

log_nH=data['log_nH']
log_Z=data['log_Z']

# plt.scatter(log_nH,log_Z,c='red',alpha=0.4)

# mask=data['O+5']==0
# data=data[mask]

# log_nH=data['log_nH']
# log_Z=data['log_Z']

# plt.scatter(log_nH,log_Z,c='green')


# def col_vs_nH(i):

#     Z=log_Z[i]
#     mask=log_Z==float(Z)
#     data1=data[mask]
#     nH=data1['log_nH']
#     log_col_den=log10(data1['O+5'])

#     plt.scatter(nH,log_col_den,label=f'{Z}',marker='D') 

def col_vs_nH(i):

    Z=log_Z[i]
    mask=log_Z==float(Z)
    data1=data[mask]
    nH=data1['log_nH']
    log_col_den=log10(data1['O+5'])

    plt.scatter(nH,log_col_den,label=f'{Z}',marker='D') 

# Z=arange(-3,2.01,0.2)
# idx=[indexOf(log_Z,round(z,2)) for z in Z]

# for i in idx:
#     col_vs_nH(i)

# plt.legend()


# plt.show()
# Z=arange(-3,2.01,0.2)
# idx=[indexOf(log_Z,round(z,2)) for z in Z]

# for i in idx:
#     col_vs_nH(i)

# plt.legend()


# plt.show()


Z_interp=unique(log_Z)
nH_interp=unique(log_nH)



def interp_func(z,kind='quadratic'):

    mask=log_Z==z
    data1=data[mask]
    nH=data1['log_nH']
    col_den=data1['O+5']

    for i,j in enumerate(col_den):
        if j==0:
            break
    
    log_nH_func=nH[:i]
    log_col_den=log10(col_den[:i])

    if 0.64 <= z <= 0.88:
        f=interp1d(log_nH_func[:-1],log_col_den[:-1],fill_value='extrapolate',kind=kind)
    
    else:
        f=interp1d(log_nH_func,log_col_den,fill_value='extrapolate',kind=kind)


    return f

# i=arange(0,len(Z_interp),1)


# for j in i:

#     plt.figure()

#     f=interp_func(Z_interp[j])
#     nH=arange(-5,0.01,0.005)
    
#     plt.scatter(nH,f(nH),label='interp')
#     col_vs_nH(indexOf(log_Z,Z_interp[j]))
#     plt.title(f'{Z_interp[j]}')
#     plt.legend()

# plt.show()


x=array(list(nH_interp)*len(Z_interp))
y=zeros(len(x))

col_den=zeros((len(Z_interp),len(nH_interp)))

for i,z in enumerate(Z_interp):

    f=interp_func(z)
    col_den[i]=f(nH_interp)
    y[len(nH_interp)*i:len(nH_interp)*(i+1)]=z


col_den=col_den.flatten()

tab=Table()

tab.add_column(x,name='log_nH')
tab.add_column(y,name='log_Z')
tab.add_column(col_den,name='O+5')

tab.write('Data/col_den_OVI.fits')

fig=plt.figure()
ax=plt.axes(projection ='3d')

ax.scatter(x,y,col_den)
ax.set_xlabel('nH')
ax.set_ylabel('Z')
ax.set_zlabel('col den')

plt.show()

