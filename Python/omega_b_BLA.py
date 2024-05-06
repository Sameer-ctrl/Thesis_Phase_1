from  numpy import *
from scipy.constants import *
from astropy.io import ascii
import matplotlib.pyplot as plt
from astropy.table import Table
from scipy.integrate import quad
from roman import toRoman
import matplotlib
import os

plt.style.use('../Python/my_style.mpl')
# matplotlib.rcParams['backend']='TkAgg'


qso_list=loadtxt('../Python/Data/qso_list.txt',dtype=str)
qso_dict=dict(zip(qso_list[:,1],qso_list[:,0]))

data=ascii.read('../Cloudy_runs/ionisation_modelling_sol.txt')
qso_all=data['qso']
z_abs_all=data['z_abs']

v_z=lambda z : 3e5*(((1+round(z,6))**2-1)/((1+round(z,6))**2+1))  # v at z
err_vz=lambda z,z_err: 4*3e5*((1+round(z,6))/(((1+round(z,6))**2)+1)**2)*round(z_err,6)
z_v=lambda v : sqrt((1+((v)/3e5))/(1-((v)/3e5)))-1      # z at v

def integrand_X(z,omega_m=0.31,omega_lambda=0.69):
    return ((1+z)**2)/sqrt(omega_lambda+(omega_m*((1+z)**3)))

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

    def __init__(self,qso,z_abs,cont_mark='*',fix_param_mark='A'):

        file=f'../VPfit/{qso}/z={z_abs:.6f}/fit_params.txt'
        
        with open(file) as f:
            text=f.read()
            text=text.replace(fix_param_mark,'')
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

        os.remove('temp_param_file.txt')

def redshift_path_qo(qso,wave_min=1220,v_lim=5000):

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

# dz=0

# qso=['3c263','pks0637', 'pg1424', 'pg0003', 'pg1216', 's135712', '1es1553', 'sbs1108', 'pg1222', 'pg1116', 'h1821', 'pg1121', 'pks0405','he0056', 'rxj0439', 'uks0242', 'pg1259', 'pks1302', '3c57', 'p1103', 'phl1811', 'pg0832']

# for q in qso:
#     dz+=redshift_path_qo(q)[1]

# print(dz)
# quit()

def redshift_path_lambda_CDM(qso,wave_min=1220,v_lim=5000):

    file_systems=open(f'Data/IGM_Danforth_Data/Systems/{qso}_igm-systems.txt','r')
    z_em=float(file_systems.readlines()[16].split(' ')[1])

    data=ascii.read(f'Data/IGM_Danforth_Data/Cont_norm_spectra/{qso}_cont_norm.asc')
    wave=data['WAVE']

    excluded_wave=ascii.read(f'Data/IGM_Danforth_Data/Excluded_wavelengths/{qso}_excluded_wavelength.asc')

    rest_wave=1215.6701

    wave_l=excluded_wave['WAVE1']
    wave_l=sort(wave_l)
    wave_r=excluded_wave['WAVE2']
    wave_r=sort(wave_r)

    z_l=(wave_l-rest_wave)/rest_wave
    z_r=(wave_r-rest_wave)/rest_wave

    z_lim=z_v(v_z(z_em)-v_lim)

    zmax=round(min([z_lim,(wave[-1]-rest_wave)/rest_wave]),3)
    zmin=round((wave_min-rest_wave)/rest_wave,6)

    delta_X=0
        
    for i in range(len(excluded_wave)):

        if z_r[i] >= zmin > z_l[i]:
            a=0
            break

        elif z_l[i]>=zmin:
            a=1
            break

    for j in range(i,len(excluded_wave)):

        if j==len(excluded_wave)-1:
            delta_X+=quad(integrand_X,z_r[j],zmax)[0]

        else:
            delta_X+=quad(integrand_X,z_r[j],z_l[j+1])[0]

    if a==1:
        delta_X+=quad(integrand_X,zmin,z_l[i])[0]

    return round(delta_X,3)


def f_H(b,x=100):

    'f_H=(HI+HII)/HI ~ HII/HI'

    logT=log10(mH*((b*1000)**2)*(1/(2*Boltzmann))*((x/100)**2))
    log_fH=(5.4*logT)-(0.33*(logT**2))-13.9   # ref. Sutherland & Dopita 1993, Richter et al. 2004, Tracing baryons in WHIM using BLAs

    return log_fH,logT

def omega_BLA(qso,b,N,err_b,err_N):

    if not isinstance(qso, (list, tuple, type(array([])))):
        qso=array([qso])

    if not isinstance(b, (list, tuple, type(array([])))):

        b=array([b])
        N=array([N])
        err_b=array([err_b])
        err_N=array([err_N])

    if isinstance(b, (list, tuple)):
        b=array(b)
        
    if isinstance(N, (list, tuple)):
        N=array(N)

    if isinstance(err_b, (list, tuple)):
        err_b=array(err_b)
        
    if isinstance(err_N, (list, tuple)):
        err_N=array(err_N)
            
    delta_X=sum([redshift_path_lambda_CDM(q) for q in unique(qso)])
    fH, logT = f_H(b)
    NH=sum(10**(fH+N))

    omega_BLA_val=A*(NH/delta_X)

    err_logT=2*(err_b/b)
    err_fH=(5.4-(0.66*(logT)))*err_logT

    err_omega_BLA=log(10)*(A/delta_X)*sqrt(sum((10**(2*(fH+N)))*(err_fH**2+err_N**2)))

    return omega_BLA_val,err_omega_BLA

b_BLA_th=40

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
            abs_system('pks0405',0.167125),
            abs_system('he0056',0.043265),
            # abs_system('pg1216',0.006328),
            abs_system('3c263',0.063397),
            abs_system('pg1222',0.054479),
            abs_system('rxj0439',0.005568),                        
            # abs_system('uks0242',0.063850),
            abs_system('pg1259',0.046284),
            abs_system('pks1302',0.094839),
            abs_system('3c57',0.077430),
            abs_system('p1103',0.003934),
            abs_system('phl1811',0.080928),
            # abs_system('pg0832',0.017505,cont_mark='^',fix_param_mark='B')
           ]

BLA_dict={}
z_abs=[]

for a in absorbers:
    BLA_obj=a.BLA_obj
    b=BLA_obj.b

    if len(b)>0:
        BLA_dict[a.qso]=[[],[],[],[]]
        z_abs.append(a.z_abs)

z_abs=unique(z_abs)

for a in absorbers:
    BLA_obj=a.BLA_obj
    b=BLA_obj.b
    N=BLA_obj.logN

    if len(b)>0:

        b_val=[]
        b_err=[]
        N_val=[]
        N_err=[]

        for i in range(len(b)):
            b_val.append(b[i][0])
            b_err.append(b[i][1])
            N_val.append(N[i][0])
            N_err.append(N[i][1])

        BLA_dict[a.qso][0]+=b_val
        BLA_dict[a.qso][1]+=N_val
        BLA_dict[a.qso][2]+=b_err
        BLA_dict[a.qso][3]+=N_err

qso=list(BLA_dict.keys())
# for q in qso:
#     print(q,redshift_path_lambda_CDM(q))

h=0.7
H0=100*h    #km/s/ Mpc
mu=1.3
mH=(proton_mass+electron_mass)   #kg

rho_c=2*(H0**2)*(1/(8*pi*gravitational_constant))
A=(mu*mH*H0)/(rho_c*speed_of_light)*(parsec*1e7)   #cm^2

omega_BLA_los=[]
omega_BLA_los_err=[]

b_all=[]
b_all_err=[]
N_all=[]
N_all_err=[]

for i,q in enumerate(qso):
    
    b_los, N_los, b_los_err, N_los_err = BLA_dict[q]

    b_all+=b_los
    N_all+=N_los
    
    b_all_err+=b_los_err
    N_all_err+=N_los_err

    val, err = omega_BLA(q,b_los,N_los,b_los_err,N_los_err)

    omega_BLA_los.append(val*100)
    omega_BLA_los_err.append(err*100)

    print(q, round(val*100,2))

# qso=['3c263', 'pks0637', 'pg1424', 'pg0003', 'pg1216', 's135712', '1es1553', 'pg1222', 'pg1116', 'h1821', 'pg1121', 'pks0405']
# b_all=[87, 162, 45, 46, 40, 63, 40, 66, 64, 52, 53, 46, 51, 64, 52, 43, 71, 63, 84, 62, 60, 56]
# N_all=[13.49, 13.6, 15.01, 14.61, 13.49, 14.2, 14.1, 13.37, 14.17, 15.1, 13.15, 15.01, 13.88, 13.54, 14.34, 15.43, 13.6, 13.68, 13.64, 13.48, 14.34, 13.09]
# b_all_err=[10, 21, 1, 4, 3, 0, 4, 10, 3, 3, 10, 4, 1, 19, 4, 1, 14, 3, 13, 11, 6, 9]
# N_all_err=[0.06, 0.06, 0.02, 0.07, 0.02, 0.02, 0.05, 0.05, 0.04, 0.05, 0.18, 0.16, 0.01, 0.11, 0.05, 0.04, 0.23, 0.02, 0.11, 0.06, 0.09, 0.06]


omega_BLA_all, omega_BLA_all_err = omega_BLA(qso,b_all,N_all,b_all_err,N_all_err)
print(f'\u03A9 = {round(omega_BLA_all*100,2)} \u00B1 {round(omega_BLA_all_err*100,2)}')

# plt.figure()

# plt.hist(omega_BLA_los,bins='auto',histtype='step')
# plt.xlabel('$\mathbf{\Omega_b(BLA) \ [\\times 10^{-2}]}$')
# plt.vlines(omega_BLA_all*100,0,8,ls='--',color='red')


# plt.figure()

# plt.hist(omega_BLA_los_err,bins='auto',histtype='step')
# plt.xlabel('$\mathbf{\Delta(\Omega_b(BLA)) \ [\\times 10^{-2}]}$')
# plt.vlines(omega_BLA_all_err*100,0,8,ls='--',color='red')

# plt.figure()

# bin_size=0.05
# bins=int((max(z_abs)-min(z_abs))/bin_size)

# plt.hist(z_abs,histtype='step',bins=bins)
# plt.xlabel('$\mathbf{z}$')


# plt.figure()
# plt.errorbar(qso,omega_BLA_los,yerr=omega_BLA_los_err,fmt='o',capsize=3,color='red')
# plt.hlines(omega_BLA_all*100,qso[0],qso[-1],lw=3,color='black')
# plt.hlines((omega_BLA_all-omega_BLA_all_err)*100,qso[0],qso[-1],ls='--',lw=2,color='black')
# plt.hlines((omega_BLA_all+omega_BLA_all_err)*100,qso[0],qso[-1],ls='--',lw=2,color='black')


# plt.figure()

# plt.hist(b_all)

# plt.figure()

# plt.hist(N_all)

# plt.show()










































# qso=['3c263', 'pks0637', 'pks0637', 'pg1424', 'pg0003', 'pg0003', 'pg0003', 'pg1216', 's135712', '1es1553', 'sbs1108', 'pg1222', 'pg1116', 'h1821', 'h1821', 'pg1121', 'pks0405']
# b=array([87.0, 162.0, 46.0, 29.0, 63.0, 40.0, 64.0, 52.0, 46.0, 51.0, 16.0, 52.0, 71.0, 63.0, 84.0, 60.0, 26.0])
# N=array([13.49, 13.6, 14.61, 15.44, 14.2, 14.1, 14.17, 15.1, 15.01, 13.88, 15.79, 14.34, 13.6, 13.68, 13.64, 14.34, 13.46])
