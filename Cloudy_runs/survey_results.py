from numpy import *
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from roman import toRoman
from io import StringIO
from astropy.io import ascii


plt.style.use('../Python/my_style.mpl')

b_BLA_th=40

qso_list=loadtxt('../Python/Data/qso_list.txt',dtype=str)
qso_dict=dict(zip(qso_list[:,1],qso_list[:,0]))

data=ascii.read('ionisation_modelling_sol.txt')

qso_all=data['qso']
z_abs_all=data['z_abs']


def ion_label(ion,ion_font_size=25,radicle_font_size=17):

    a=ion.split('+')

    if len(a)>1:

        if a[1]!='':
            return f'{{\\fontsize{{{ion_font_size}pt}}{{3em}}\selectfont{{}}$\mathbf{{{a[0]}}}$}} {{\\fontsize{{{radicle_font_size}pt}}{{3em}}\selectfont{{}}$\mathbf{{{toRoman(int(a[1])+1)}}}$}}'

        else:
            return f'{{\\fontsize{{{ion_font_size}pt}}{{3em}}\selectfont{{}}$\mathbf{{{a[0]}}}$}} {{\\fontsize{{{radicle_font_size}pt}}{{3em}}\selectfont{{}}$\mathbf{{{toRoman(2)}}}$}}'

    else:

        return f'{{\\fontsize{{{ion_font_size}pt}}{{3em}}\selectfont{{}}$\mathbf{{{a[0]}}}$}} {{\\fontsize{{{radicle_font_size}pt}}{{3em}}\selectfont{{}}$\mathbf{{{toRoman(1)}}}$}}'

v_z=lambda z : 3e5*(((1+round(z,6))**2-1)/((1+round(z,6))**2+1))
err_vz=lambda z,z_err: 4*3e5*((1+round(z,6))/(((1+round(z,6))**2)+1)**2)*round(z_err,6)


class ion():

    def __init__(self,name,v,b,logN):

        self.ion=name
        self.v=[x for x in v]
        self.b=[x for x in b]
        self.logN=[x for x in logN]
        self.comp=len(v)


class abs_system():

    def __init__(self,qso,z_abs,cont_mark='*'):

        file=f'../VPfit/{qso}/z={z_abs}/fit_params.txt'

        with open(file) as f:
            text=f.read()
            text=text.replace('A','')
            # print(text)
            # print('\n')

        with open('temp_param_file.txt','w') as f:
            f.write(text)

        param_file=loadtxt('temp_param_file.txt',dtype=str,comments=('>','#'))

        ions=param_file[:,0]
        mask=[]

        for i in ions:
            mask.append(cont_mark not in i)  

        ions=ions[mask]
        z=param_file[:,1].astype(float)[mask]
        z_err=param_file[:,2].astype(float)[mask]
        b=param_file[:,3].astype(float)[mask]
        b_err=param_file[:,4].astype(float)[mask]
        logN=param_file[:,5].astype(float)[mask]
        logN_err=param_file[:,6].astype(float)[mask]

        v_abs=v_z(z_abs)
        v=zeros(len(ions))
        v_err=zeros(len(ions))

        for i in range(len(ions)):

            v[i]=round(v_z(z[i])-v_abs)
            v_err[i]=round(err_vz(z[i],z_err[i]))

        ions_all=unique(ions)
        ion_obj_dict={}

        for i in ions_all:
            mask=ions==i
            v_ion=v[mask]
            v_err_ion=v_err[mask]
            b_ion=b[mask]
            b_err_ion=b_err[mask]
            logN_ion=logN[mask]
            logN_err_ion=logN_err[mask]

            v_obj=[]
            b_obj=[]
            logN_obj=[]

            for j in range(len(v_ion)):
                v_obj.append([v_ion[j],v_err_ion[j]])
                b_obj.append([round(b_ion[j]),round(b_err_ion[j])])
                logN_obj.append([round(logN_ion[j],2),round(logN_err_ion[j],2)])

            obj=ion(i,v=v_obj,b=b_obj,logN=logN_obj)
            ion_obj_dict[i]=obj

        v_BLA_obj=[]
        b_BLA_obj=[]
        logN_BLA_obj=[]

        v_HI=ion_obj_dict['HI'].v
        b_HI=ion_obj_dict['HI'].b
        logN_HI=ion_obj_dict['HI'].logN

        for i,b_val in enumerate(b_HI):

            if b_val[0]>=b_BLA_th:
                v_BLA_obj.append(v_HI[i])
                b_BLA_obj.append(b_HI[i])
                logN_BLA_obj.append(logN_HI[i])
    
        self.BLA_obj=ion('HI',v=v_BLA_obj,b=b_BLA_obj,logN=logN_BLA_obj)

        self.qso_label=qso_dict[qso]
        self.z_abs=z_abs
        self.ion_obj=ion_obj_dict
        self.ions=ions_all[ions_all!='HI']
        self.n_ions=len(self.ions)


        mask=logical_and(qso_all==qso,z_abs_all==z_abs)
        data_abs=data[mask]

        NHi_abs=data_abs['NHi']

        sol={}

        for n in unique(NHi_abs):

            NHi_mask=NHi_abs==n

            nH_abs_masked=data_abs['nH'][NHi_mask]
            err_nH_abs_masked=data_abs['err_nH'][NHi_mask]
            Z_abs_masked=data_abs['Z'][NHi_mask]
            err_Z_abs_masked=data_abs['err_Z'][NHi_mask]
            case_abs_masked=data_abs['case'][NHi_mask]

            for i,c in enumerate(case_abs_masked):
                
                if c=='exc':
                    sol_exc=[nH_abs_masked[i],err_nH_abs_masked[i],Z_abs_masked[i],err_Z_abs_masked[i]]

                else:
                    sol_inc=[nH_abs_masked[i],err_nH_abs_masked[i],Z_abs_masked[i],err_Z_abs_masked[i]]

            sol[n]=[sol_exc,sol_inc]
    
        self.ion_modelling_sol=sol


absorbers=[
            abs_system('3c263',0.140756),
            abs_system('pks0637',0.161064),
            abs_system('pks0637',0.417539),
            abs_system('pg1424',0.147104),
            abs_system('pg0003',0.347579),                        
            abs_system('pg0003',0.386089),
            abs_system('pg0003',0.421923),
            abs_system('pg1216',0.282286),
            abs_system('s135712',0.097869),
            abs_system('1es1553',0.187764),
            abs_system('sbs1108',0.463207),
            abs_system('pg1222',0.378389),
            abs_system('pg1116',0.138527),
            abs_system('h1821',0.170006),
            abs_system('h1821',0.224981),
            abs_system('pg1121',0.192393),
            abs_system('pks0405',0.167125)
           ]

# for a in absorbers:

#     BLA=a.BLA_obj
#     Ovi=a.ion_obj['OVI']
    
#     v_BLA=[x[0] for x in BLA.v]
#     v_err_BLA=[x[1] for x in BLA.v]

#     v_Ovi=[x[0] for x in Ovi.v]
#     v_err_Ovi=[x[1] for x in Ovi.v]

#     plt.title(f'{a.qso_label} (z_abs={a.z_abs})')
#     plt.errorbar(v_BLA,10*ones(len(v_BLA)),0,v_err_BLA,label='BLA',fmt='o',capsize=3)
#     plt.errorbar(v_Ovi,9*ones(len(v_Ovi)),0,v_err_Ovi,label='OVI',fmt='o',capsize=3)
#     plt.legend()
#     plt.ylim(0,11)
#     plt.show()

ion_model_sol=[a.ion_modelling_sol for a in absorbers]

NHi=zeros(int(len(data)/2))

nH_exc=zeros(int(len(data)/2))
err_nH_exc=zeros(int(len(data)/2))
Z_exc=zeros(int(len(data)/2))
err_Z_exc=zeros(int(len(data)/2))

nH_inc=zeros(int(len(data)/2))
err_nH_inc=zeros(int(len(data)/2))
Z_inc=zeros(int(len(data)/2))
err_Z_inc=zeros(int(len(data)/2))

k=0

for i,abs_sol in enumerate(ion_model_sol):

    for n in abs_sol:
        exc,inc=abs_sol[n]
        
        NHi[k]=n

        nH_exc[k],err_nH_exc[k],Z_exc[k],err_Z_exc[k]=exc
        nH_inc[k],err_nH_inc[k],Z_inc[k],err_Z_inc[k]=inc

        k+=1

# plt.figure(figsize=(16,10))

# plt.subplot(121)
# plt.title(f'$\mathbf{{Excluding \ }}$'+ion_label('O+5',ion_font_size=20,radicle_font_size=15))
# # plt.scatter(nH_exc,Z_exc)
# plt.errorbar(nH_exc,Z_exc,xerr=err_nH_exc,yerr=err_Z_exc,fmt='o',capsize=3)
# plt.xlabel(r'$\mathbf{log \ n_H}$')
# plt.ylabel(r'$\mathbf{log \ Z}$')
# # plt.colorbar()

# plt.subplot(122)
# plt.title(f'$\mathbf{{Including \ }}$'+ion_label('O+5',ion_font_size=20,radicle_font_size=15))
# # plt.scatter(nH_inc,Z_inc)
# plt.errorbar(nH_inc,Z_inc,xerr=err_nH_inc,yerr=err_Z_inc,fmt='o',capsize=3)
# plt.xlabel(r'$\mathbf{log \ n_H}$')
# plt.ylabel(r'$\mathbf{log \ Z}$')
# # plt.colorbar()

# plt.subplots_adjust(wspace=0.34)
# plt.savefig('Z_vs_nH.png')


# plt.figure(figsize=(16,10))

# plt.subplot(121)
# plt.title(f'$\mathbf{{Excluding \ }}$'+ion_label('O+5',ion_font_size=20,radicle_font_size=15))
# # plt.scatter(nH_exc,NHi)
# plt.errorbar(nH_exc,NHi,xerr=err_nH_exc,fmt='o',capsize=3)
# plt.xlabel(r'$\mathbf{log \ n_H}$')
# plt.ylabel(r'$\mathbf{log \ N(Hi)}$')

# plt.subplot(122)
# plt.title(f'$\mathbf{{Including \ }}$'+ion_label('O+5',ion_font_size=20,radicle_font_size=15))
# # plt.scatter(nH_inc,NHi)
# plt.errorbar(nH_inc,NHi,xerr=err_nH_inc,fmt='o',capsize=3)
# plt.xlabel(r'$\mathbf{log \ n_H}$')
# plt.ylabel(r'$\mathbf{log \ N(Hi)}$')

# plt.subplots_adjust(wspace=0.34)
# plt.savefig('NHi_vs_nH.png')


plt.figure(figsize=(16,10))

plt.subplot(121)
plt.title(f'$\mathbf{{Excluding \ }}$'+ion_label('O+5',ion_font_size=20,radicle_font_size=15))
plt.errorbar(Z_exc,NHi,xerr=err_Z_exc,fmt='o',capsize=3)
plt.xlabel(r'$\mathbf{log \ Z}$')
plt.ylabel(r'$\mathbf{log \ N(Hi)}$')

plt.subplot(122)
plt.title(f'$\mathbf{{Including \ }}$'+ion_label('O+5',ion_font_size=20,radicle_font_size=15))
plt.errorbar(Z_inc,NHi,xerr=err_Z_inc,fmt='o',capsize=3)
plt.xlabel(r'$\mathbf{log \ Z}$')
plt.ylabel(r'$\mathbf{log \ N(Hi)}$')

plt.subplots_adjust(wspace=0.34)
plt.savefig('NHi_vs_Z.png')

# plt.figure()

# plt.subplot(121)
# plt.hist(nH_exc,bins=4)
# plt.xlabel(r'$\mathbf{log \ n_H}$')

# plt.subplot(122)
# plt.hist(nH_inc,bins=4)
# plt.xlabel(r'$\mathbf{log \ n_H}$')

plt.show()























































quit()

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
        self.b_OVI_err=b_OVI_err
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
b_H_err=zeros(len(absorbers))
n_H=zeros(len(absorbers))
n_H_err=zeros(len(absorbers))
b_OVI=zeros(len(absorbers))
b_OVI_err=zeros(len(absorbers))
n_OVI=zeros(len(absorbers))
n_ions=zeros(len(absorbers))
is_BLA=zeros(len(absorbers))

for i,a in enumerate(absorbers):
    b_H[i]=a.b_H
    b_H_err[i]=a.b_H_err
    b_OVI[i]=a.b_OVI
    b_OVI_err[i]=a.b_OVI_err
    n_H[i]=a.nH
    n_H_err[i]=a.nH_err
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


# plt.figure(figsize=(8,5))

# plt.scatter(b_H,b_OVI,c='red')
# plt.errorbar(b_H,b_OVI,yerr=b_OVI_err,xerr=b_H_err,c='red',fmt='o',capsize=3)
# plt.plot(array([0,150]),array([0,150]),ls='--')
# plt.plot(array([0,165]),array([0,165/4]),ls='--')
# plt.vlines(40,-10,170,ls='--',color='green')
# plt.xlabel(f'$\mathbf{{b(}}$'+ion_label('H',ion_font_size=20,radicle_font_size=13)+r'$\mathbf{) \ [ \ km \ {s}^{-1}}]$',labelpad=10,fontsize=20)
# plt.ylabel(f'$\mathbf{{b(}}$'+ion_label('O+5',ion_font_size=20,radicle_font_size=13)+r'$\mathbf{) \ [ \ km \ {s}^{-1}}]$',labelpad=10,fontsize=20)
# plt.ylim(bottom=0,top=160)
# plt.xlim(0,170)
# plt.savefig('bHi_vs_BOvi_1.png')
# plt.show()
# quit()


# plt.figure(figsize=(13,7))

# plt.subplot(121)
# plt.title(r'$\mathbf{Column \ density}$',fontsize=25)
# plt.hist(n_H,bins=7)
# plt.xlabel(f'$\mathbf{{log \ [N(}}$'+ion_label('H')+r'$\mathbf{) \ {cm}^{-2}}]$',labelpad=15)
# plt.ylabel(r'$\mathbf{No. \ of \ absorbers}$',labelpad=15)

# plt.subplot(122)
# plt.title(r'$\mathbf{Doppler \ parameter}$',fontsize=25)
# plt.hist(b_H,bins='auto')
# plt.xlabel(f'$\mathbf{{b \ (}}$'+ion_label('H')+r'$\mathbf{) \ [km \ {s}^{-1}}]$',labelpad=15)
# plt.ylabel(r'$\mathbf{No. \ of \ absorbers}$',labelpad=15)


# plt.subplots_adjust(wspace=0.34)

# plt.savefig('NHi_vs_bHi.png')


plt.figure(figsize=(8,5))

plt.hlines(40,13,16,ls='--',color='black',lw=3)
# plt.scatter(n_H,b_H,label='H',s=70)
plt.errorbar(n_H,b_H,xerr=n_H_err,yerr=b_H_err,fmt='o',capsize=3,c='red')
# plt.plot([42,100],f(array([42,100]),fit[0][0],fit[0][1]),ls='--')

# plt.xlabel(r'$\mathbf{b} \ \mathbf{(km \ \ s^{-1})}$',labelpad=15)
# plt.ylabel(f'$\mathbf{{log \ [N(}}$'+ion_label('H')+r'$\mathbf{) \ {cm}^{-2}}]$',labelpad=15)
plt.xlim(13.35,15.9)

# plt.scatter(n_H,b_H,label='H',c=is_BLA)
# plt.plot(f(array([42,100]),fit[0][0],fit[0][1]),[42,100],ls='--')
# plt.ylabel(r'$\mathbf{b} \ \mathbf{(km \ \ s^{-1})}$',labelpad=15)
plt.ylabel(f'$\mathbf{{b \ (}}$'+ion_label('H')+r'$\mathbf{) \ [km \ {s}^{-1}}]$',labelpad=15)
plt.xlabel(f'$\mathbf{{log \ [N(}}$'+ion_label('H')+r'$\mathbf{) \ {cm}^{-2}}]$',labelpad=15)
plt.savefig('NHi_vs_bHi.png')
plt.show()



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


        