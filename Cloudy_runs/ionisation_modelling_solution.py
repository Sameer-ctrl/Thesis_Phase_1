import os
from numpy import *
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from roman import toRoman
from io import StringIO
from astropy.io import ascii
from time import sleep



plt.style.use('../Python/my_style.mpl')


cwd=os.getcwd()
cwd_split=cwd.split('/')

cloudy_path_desktop='/home/sameer/Documents/Sameer/c22.02/source/cloudy.exe'
cloudy_path_workstation='/home/sameer/Sameer/c22.02/source/cloudy.exe'
cloudy_path_pc='/home/sameer/cloudy/source/sys_gcc/cloudy.exe'

if cwd_split[3]=='Thesis_Phase_1':
    cloudy_path=cloudy_path_pc

elif cwd_split[3]=='Sameer':
    cloudy_path=cloudy_path_workstation

else:
    cloudy_path=cloudy_path_desktop


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
        self.qso=qso
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


def cloudy_run_sol():

    for a in absorbers:

        qso_los=a.qso
        z_abs=a.z_abs
        ion_model_sol=a.ion_modelling_sol

        for N in ion_model_sol:

            nH=ion_model_sol[N][0][0]
            Z=ion_model_sol[N][0][2]
            run_name=f'{qso_los}_{z_abs}_{N}'
            
            with open('template.in') as f:
                data=f.read()
                data=data.replace('z_abs',f'{z_abs}')
                data=data.replace('nH',f'{nH}')
                data=data.replace('Z',f'{Z}')
                data=data.replace('NHi',f'{N}')
                data=data.replace('run_name',f'{run_name}')

            with open(f'ionisation_modelling_solution/{run_name}.in','w') as f:
                f.write(data)
            
            os.chdir(f'ionisation_modelling_solution')
            
            cloudy_run_command=f'{cloudy_path} -r {run_name}'

            print(f'\nCloudy running in ionisation_modelling_solution : {qso_los} : {z_abs} : {N} \n')
            os.system(cloudy_run_command)

            sleep(1)

            os.remove(f'{run_name}.out')

            print('----------- Writing output files ------------- \n')


            with open(f'{run_name}_col_density.txt') as f:
                data=f.read()
                data=data.replace('#column density ','')

            with open(f'{run_name}_col_density.txt','w') as f:
                f.write(data)

            col_density=ascii.read(f'{run_name}_col_density.txt')

            temp_file=f'{run_name}_temp.txt'

            data_temp=genfromtxt(temp_file,delimiter=[11,11,9,10,10])

            Te=data_temp[:,1]
            d2t_dr2=data_temp[:,4]

            log_Te=zeros(len(col_density))
            k=0

            for i,j in enumerate(d2t_dr2):
                if j==0:
                    log_Te[k]=round(log10(Te[i]),3)
                    k+=1

            col_density.add_column(log_Te,name='log_Te')
            col_density.add_column([nH],name='log_nH')
            col_density.add_column([Z],name='log_Z')
            col_density.add_column([qso_los],name='qso')
            col_density.add_column([z_abs],name='z_abs')
            col_density.add_column([N],name='NHi')

            os.remove(temp_file)

            col_density.write(f'output_fits/{run_name}_col_density_param.fits',overwrite=True)
            ascii.write(col_density,f'output_txt/{run_name}_col_density_param.txt',format='ecsv',overwrite=True)
                

            print('----------- Output files written -------------')
            os.chdir(cwd)


cloudy_run_sol()
