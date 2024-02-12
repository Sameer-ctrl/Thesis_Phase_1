from numpy import *
import matplotlib.pyplot as plt


class BLA():

    def __init__(self,qso,z_abs,nH,b_Lya,n_OVI,b_OVI,n_ions,bla):

        self.los=qso
        self.z=z_abs
        self.name=f'{qso}_z={z_abs}'
        self.nH=nH
        self.b_Lya=b_Lya
        self.n_OVI=n_OVI
        self.b_OVI=b_OVI
        self.n_ions=n_ions
        self.is_BLA=bla


absorbers=[
    BLA('3c263',0.140756,13.49,87,13.63,26,3,True),
    BLA('pks0637',0.161064,13.60,162,14.02,48,3,True),
    BLA('pks0637',0.417539,14.61,46,14.19,42,3,True),
    BLA('pg1424',0.147104,15.44,29,13.73,16,4,False),
    BLA('pg0003',0.347579,14.20,63,14.25,30,5,True),                        
    BLA('pg0003',0.386089,14.1,40,13.71,25,4,False),
    BLA('pg0003',0.421923,14.17,64,14.27,27,3,True),
    BLA('pg1216',0.282286,15.1,52,13.93,58,3,True),
    BLA('s135712',0.097869,15.01,46,14.30,43,5,True),
    BLA('1es1553',0.187764,13.88,51,14.23,3,3,True),
    BLA('sbs1108',0.463207,15.79,16,13.71,45,7,False),
    BLA('pg1222',0.378389,14.34,52,13.68,34,4,True),
    BLA('pg1116',0.138527,13.6,71,13.84,35,9,True),
    BLA('h1821',0.170006,13.68,63,13.94,152,3,True),
    BLA('h1821',0.224981,13.64,84,14.24,45,3,True),
    BLA('pg1121',0.192393,14.34,60,12.84,11,6,True),
    BLA('pks0405',0.167125,13.46,26,14.05,41,10,False)]

b_Lya=zeros(len(absorbers))
n_H=zeros(len(absorbers))
b_OVI=zeros(len(absorbers))
n_OVI=zeros(len(absorbers))
n_ions=zeros(len(absorbers))
is_BLA=zeros(len(absorbers))

for i,a in enumerate(absorbers):
    b_Lya[i]=a.b_Lya
    b_OVI[i]=a.b_OVI
    n_H[i]=a.nH
    n_OVI[i]=a.n_OVI
    n_ions=a.n_ions
    is_BLA[i]=a.is_BLA

plt.subplot(121)
plt.title('Lya')
plt.scatter(b_Lya,n_H,label='Lya',c=is_BLA)
plt.colorbar()

plt.subplot(122)
plt.title('OVI')
plt.scatter(b_OVI,n_OVI,label='OVI',c=is_BLA)
plt.colorbar()

# plt.legend()
plt.show()


        