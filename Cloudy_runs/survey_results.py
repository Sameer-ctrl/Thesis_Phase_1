from numpy import *
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from roman import toRoman
from io import StringIO

plt.style.use('../Python/my_style.mpl')

def ion_label(ion,ion_font_size=25,radicle_font_size=17):

    a=ion.split('+')

    if len(a)>1:

        if a[1]!='':
            return f'{{\\fontsize{{{ion_font_size}pt}}{{3em}}\selectfont{{}}$\mathbf{{{a[0]}}}$}} {{\\fontsize{{{radicle_font_size}pt}}{{3em}}\selectfont{{}}$\mathbf{{{toRoman(int(a[1])+1)}}}$}}'

        else:
            return f'{{\\fontsize{{{ion_font_size}pt}}{{3em}}\selectfont{{}}$\mathbf{{{a[0]}}}$}} {{\\fontsize{{{radicle_font_size}pt}}{{3em}}\selectfont{{}}$\mathbf{{{toRoman(2)}}}$}}'

    else:

        return f'{{\\fontsize{{{ion_font_size}pt}}{{3em}}\selectfont{{}}$\mathbf{{{a[0]}}}$}} {{\\fontsize{{{radicle_font_size}pt}}{{3em}}\selectfont{{}}$\mathbf{{{toRoman(1)}}}$}}'




class BLA():

    def __init__(self,qso,z_abs,nH,nH_err,b_H,b_H_err,n_OVI,n_OVI_err,b_OVI,b_OVI_err,n_ions,bla):

        self.los=qso
        self.z=z_abs
        self.name=f'{qso}_z={z_abs}'
        self.nH=nH
        self.nH_err=nH_err
        self.b_H=b_H
        self.b_H_err=b_H_err
        self.n_OVI=n_OVI
        self.b_OVI=b_OVI
        self.n_ions=n_ions
        self.is_BLA=bla


absorbers=[
    BLA('3c263',0.140756,13.49,0.06,87,10,13.63,0.04,26,4,3,True),
    BLA('pks0637',0.161064,13.60,0.06,162,21,14.02,0.03,48,5,3,True),
    BLA('pks0637',0.417539,14.61,0.07,46,4,14.19,0.05,42,6,3,True),
    BLA('pg1424',0.147104,15.44,0.14,29,2,13.73,0.11,16,6,4,False),
    BLA('pg0003',0.347579,14.20,0.02,63,1,14.25,0.02,30,2,5,True),                        
    BLA('pg0003',0.386089,14.1,0.05,40,4,13.71,0.06,25,4,4,False),
    BLA('pg0003',0.421923,14.17,0.04,64,3,14.27,0.02,27,1,3,True),
    BLA('pg1216',0.282286,15.1,0.05,52,3,13.93,0.05,58,9,3,True),
    BLA('s135712',0.097869,15.01,0.16,46,4,14.30,0.11,43,16,5,True),
    BLA('1es1553',0.187764,13.88,0.01,51,1,14.23,0.33,3,1,3,True),
    BLA('sbs1108',0.463207,15.79,0.11,16,1,13.71,0.07,45,10,7,False),
    BLA('pg1222',0.378389,14.34,0.05,52,4,13.68,0.24,34,13,4,True),
    BLA('pg1116',0.138527,13.6,0.23,71,14,13.84,0.02,35,3,9,True),
    BLA('h1821',0.170006,13.68,0.02,63,3,13.94,0.06,152,20,3,True),
    BLA('h1821',0.224981,13.64,0.11,84,13,14.24,0.01,45,1,3,True),
    BLA('pg1121',0.192393,14.34,0.09,60,6,12.84,0.19,11,16,6,True),
    BLA('pks0405',0.167125,13.46,0.04,26,3,14.05,0.1,41,3,10,False)]

b_H=zeros(len(absorbers))
n_H=zeros(len(absorbers))
b_OVI=zeros(len(absorbers))
n_OVI=zeros(len(absorbers))
n_ions=zeros(len(absorbers))
is_BLA=zeros(len(absorbers))

for i,a in enumerate(absorbers):
    b_H[i]=a.b_H
    b_OVI[i]=a.b_OVI
    n_H[i]=a.nH
    n_OVI[i]=a.n_OVI
    n_ions=a.n_ions
    is_BLA[i]=a.is_BLA

b_H_fit=[]
n_H_fit=[]

for i in range(len(n_H)):

    if 100>=b_H[i]>45:
        b_H_fit.append(b_H[i])
        n_H_fit.append(n_H[i])

def f(x,m,c):
    return (m*x)+c

fit=curve_fit(f,b_H_fit,n_H_fit)
print(fit[0])

plt.subplot(121)
plt.title('N(HI)')
plt.hist(n_H,bins=7)

plt.subplot(122)
plt.title('b')
plt.hist(b_H,bins='auto')
plt.show()


# plt.scatter(b_H,n_H,label='H',c=is_BLA)
# plt.plot([42,100],f(array([42,100]),fit[0][0],fit[0][1]),ls='--')
# plt.xlabel(r'$\mathbf{b} \ \mathbf{(km \ \ s^{-1})}$',labelpad=15)
# plt.ylabel(f'$\mathbf{{log \ [N(}}$'+ion_label('H')+r'$\mathbf{) \ {cm}^{-2}}]$',labelpad=15)

# plt.scatter(n_H,b_H,label='H',c=is_BLA)
# plt.plot(f(array([42,100]),fit[0][0],fit[0][1]),[42,100],ls='--')
# plt.ylabel(r'$\mathbf{b} \ \mathbf{(km \ \ s^{-1})}$',labelpad=15)
# plt.xlabel(f'$\mathbf{{log \ [N(}}$'+ion_label('H')+r'$\mathbf{) \ {cm}^{-2}}]$',labelpad=15)
# plt.show()



# plt.subplot(121)
# plt.title('H')
# plt.scatter(b_H,n_H,label='H',c=is_BLA)
# plt.colorbar()

# plt.subplot(122)
# plt.title('OVI')
# plt.scatter(b_OVI,n_OVI,label='OVI',c=is_BLA)
# plt.colorbar()

# plt.legend()
# plt.show()


        