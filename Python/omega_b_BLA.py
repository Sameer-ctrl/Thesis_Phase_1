from  numpy import *
from scipy.constants import *
import os
from astropy.io import ascii,fits
import matplotlib.pyplot as plt
from astropy.table import Table,Column

mH=(proton_mass+electron_mass)

v_z=lambda z : 3e5*(((1+round(z,6))**2-1)/((1+round(z,6))**2+1))  # v at z
err_vz=lambda z,z_err: 4*3e5*((1+round(z,6))/(((1+round(z,6))**2)+1)**2)*round(z_err,6)
z_v=lambda v : sqrt((1+((v)/3e5))/(1-((v)/3e5)))-1      # z at v

b_BLA_th=40

qso_list=loadtxt('../Python/Data/qso_list.txt',dtype=str)
qso_dict=dict(zip(qso_list[:,1],qso_list[:,0]))

data=ascii.read('../Cloudy_runs/ionisation_modelling_sol.txt')

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
        self.qso=qso
        self.z_abs=z_abs
        self.ion_obj=ion_obj_dict
        self.ions=ions_all[ions_all!='HI']
        self.n_ions=len(self.ions)


        mask=logical_and(qso_all==qso,z_abs_all==z_abs)
        data_abs=data[mask]

        NHi_abs=data_abs['NHi']
        # err_NHi_abs=data_abs['err_NHi']

        sol={}

        for n in unique(NHi_abs):

            NHi_mask=NHi_abs==n

            err_NHi_abs_masked=data_abs['err_NHi'][NHi_mask]
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

            sol[n]=[sol_exc,sol_inc,err_NHi_abs_masked[0]]
    
        self.ion_modelling_sol=sol

def redshift_path(qso,wave_min=1220,v_lim=5000):

    file_systems=open(f'Data/IGM_Danforth_Data/Systems/{qso}_igm-systems.txt','r')
    z_em=float(file_systems.readlines()[16].split(' ')[1])

    data=ascii.read(f'Data/IGM_Danforth_Data/Cont_norm_spectra/{qso}_cont_norm.asc')
    wave=data['WAVE']

    excluded_wave=ascii.read(f'Data/IGM_Danforth_Data/Excluded_wavelengths/{qso}_excluded_wavelength.asc')

    wave_l=excluded_wave['WAVE1']
    wave_r=excluded_wave['WAVE2']

    z_lim=z_v(v_z(z_em)-v_lim)

    rest_wave=1215.6701
    wave_max=min([(1+z_lim)*rest_wave,wave[-1]])

    dzb=0

    for i in range(len(excluded_wave)):

        if wave_l[i] >= wave_min:
            dzb+=(wave_r[i]-wave_l[i])/rest_wave
        
        elif wave_r[i] >= wave_min > wave_l[i]:
            dzb+=(wave_r[i]-wave_min)/rest_wave

    dzb=round(dzb,3)

    zmax=round((wave_max-rest_wave)/rest_wave,3)
    zmin=round((wave_min-rest_wave)/rest_wave,3)

    dz=zmax-zmin
    dz_unblocked=round(dz-dzb,3)

    delta_X=round(0.5*(((1+zmax-(dzb/2)))**2-((1+zmin+(dzb/2)))**2),3)
    # print(delta_X,dz_unblocked,dzb,dz)

    return delta_X,dz_unblocked,dzb,dz

def f_H(b):

    'f_H=(HI+HII)/HI ~ HII/HI'

    logT=log10(mH*((b*1000)**2)*(1/(2*Boltzmann)))
    log_fH=(5.4*logT)-(0.33*(logT**2))-13.9   # ref. Sutherland & Dopita 1993, Richter et al. 2004, Tracing baryons in WHIM using BLAs

    return log_fH


absorbers=[
            abs_system('3c263',0.140756),
            abs_system('pks0637',0.161064),
            abs_system('pks0637',0.417539),
            abs_system('pg1424',0.147104),
            abs_system('pg0003',0.347586),                        
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

qso=unique([a.qso for a in absorbers])

BLA_dict={q:[[],[]] for q in qso}


for a in absorbers:
    BLA_obj=a.BLA_obj
    b=BLA_obj.b
    logN=BLA_obj.logN

    for i in range(len(b)):
        BLA_dict[a.qso][0].append(b[i][0])
        BLA_dict[a.qso][1].append(logN[i][0])

BLA_dict.pop('sbs1108')

h=0.7
H0=100*h    #km/s/ Mpc
mu=1.3

rho_c=2*(H0**2)*(1/(8*pi*gravitational_constant))

A=(mu*mH*H0)/(rho_c*speed_of_light)*(parsec*1e7)   # cm^2

def omega_BLA(qso,b,N):

    if not isinstance(qso, (list, tuple, type(array([])))):
        qso=array([qso])

    if not isinstance(b, (list, tuple, type(array([])))):

        b=array([b])
        N=array([N])

    if isinstance(b, (list, tuple)):
        b=array(b)
        
    if isinstance(N, (list, tuple)):
        N=array(N)
            
    print(type(b),type(N))
    delta_X=sum([redshift_path(q)[0] for q in unique(qso)])
    fH=f_H(b)
    NH=sum(10**(fH+N))

    return A*(NH/delta_X)


# qso=['3c263', 'pks0637', 'pks0637', 'pg1424', 'pg0003', 'pg0003', 'pg0003', 'pg1216', 's135712', '1es1553', 'sbs1108', 'pg1222', 'pg1116', 'h1821', 'h1821', 'pg1121', 'pks0405']
# b=array([87.0, 162.0, 46.0, 29.0, 63.0, 40.0, 64.0, 52.0, 46.0, 51.0, 16.0, 52.0, 71.0, 63.0, 84.0, 60.0, 26.0])
# N=array([13.49, 13.6, 14.61, 15.44, 14.2, 14.1, 14.17, 15.1, 15.01, 13.88, 15.79, 14.34, 13.6, 13.68, 13.64, 14.34, 13.46])


omega_BLA_all=omega_BLA(qso,b,N)*100
print(omega_BLA_all)
print(mH*((40*1000)**2)*(1/(2*Boltzmann)))

omega_BLA_los=[]

for i in range(len(qso)):
    
    los=[qso[i]]
    b_los=array([b[i]])
    N_los=array([N[i]])

    omega_BLA_los.append(omega_BLA(los,b_los,N_los)*100)

plt.figure()

plt.hist(omega_BLA_los,bins='auto')
plt.vlines(omega_BLA_all,0,7,ls='--',color='red')

plt.figure()
plt.scatter(qso,omega_BLA_los)
plt.hlines(omega_BLA_all,qso[0],qso[-1])

plt.figure()

plt.plot(linspace(0,150,1000),f_H(linspace(0,150,1000)))

plt.show()

qso='pg0003'
