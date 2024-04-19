import os
from numpy import *
from astropy.io import ascii
from time import sleep
from roman import toRoman

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


def run_cloudy(run_name, hden, metal, temp, redshift, stop_nH, ions, qso, z_abs, stop_T=1000, save_temp=False, save_Hyd=False, uvb_scale = 1, uvb_Q = 18, He_abun=0.08156498, delete_out_file=False, delete_temp_file=True, miscalleneous_command=''):

    grid_parameters=[]

    if len(hden)>1:
        hden_line = f'hden -5 log vary \ngrid range from {hden[0]} to {hden[1]} with {hden[2]} dex step \n'
        grid_parameters.append('log_nH')

    else:
        hden_line = f'hden {hden[0]} log \n'

    if len(metal)>1:
        metal_line = f'metals -1 log vary \ngrid range from {metal[0]} to {metal[1]} with {metal[2]} dex step \n'
        grid_parameters.append('log_Z')

    else:
        metal_line = f'metals {metal[0]} log \n'

    if temp!=None:

        if len(temp)>1:
            temp_line = f'constant temperature 5 log vary \ngrid range from {temp[0]} to {temp[1]} with {temp[2]} dex step \n'
            grid_parameters.append('log_T')

        else:
            temp_line = f'constant temperature {temp[0]} log \n'

    else:
        temp_line=''

    n=len(grid_parameters)

    if n>=1:
        grid_file=f'{run_name}_grid.txt'
        save_grid=f'save grid "{grid_file}" last no hash \n'

        # if temp==None:
        #     n_models=((len(range(*hden)))+1)*((len(range(*metal)))+1)
        
        # else:
        #     n_models=((len(range(*hden)))+1)*((len(range(*metal)))+1)*((len(range(*temp)))+1)
    
    else: 
        save_grid=f''


    if save_Hyd==True:
        hyd_cond=f'{run_name}_hydrogen.txt'
        save_hyd=f'save hydrogen conditions "{hyd_cond}" last no clobber \n'

    else:
        save_hyd=''


    if save_temp==True:
        temp_file=f'{run_name}_temp.txt'
        save_Temp=f'save temperature "{temp_file}" no hash \n'

    else:
        save_Temp=''
    
    col_density=f'{run_name}_col_density.txt'
    temp_file=f'{run_name}_temp.txt'


    'UV background model : KS19'

    uv_b=f'TABLE KS19 redshift = {redshift} scale = {uvb_scale} Q = {uvb_Q} \n'
    abundance=f'abundances "solar_GASS10.abn" \nelement helium abundance {He_abun} linear \n'
    stop_criteria_nH=f'stop column density {stop_nH} neutral H \nstop temperature {stop_T} K linear \n'
    save_col_den=f'save species column density "{col_density}" no hash \n'

    for i in ions:
        save_col_den+=f'"{i}" \n'

    lines=[uv_b,abundance,hden_line,metal_line,temp_line,stop_criteria_nH,miscalleneous_command,save_grid,save_Temp,save_hyd,save_col_den]

    if len(metal==1):
        path=f'{qso}/z={z_abs}/logZ={metal[0]}/{run_name}'
    
    else:
        path=f'{qso}/z={z_abs}/{run_name}'

    os.mkdir(path)
    file=open(f'{path}/{run_name}.in','w+')
    file.writelines(lines)
    file.write('end')
    file.close()


    file=open(f'{path}/{run_name}.in','r')

    print('-------------Cloudy input file------------- \n')
    print(file.read())
    print('\n--------------------------------------- \n')

    if n>=1:
        print(f'Grid parameters : {grid_parameters}')
        # print(f'No. of models : {n_models} \n')

    # a=input('Check the input file for cloudy. \nPress Y to continue.\n')

    if True:
    # if a=='y':

        os.chdir(f'{path}')

        cloudy_run_command=f'{cloudy_path} -r {run_name}'

        print(f'\nCloudy running in {path} ... \n')
        os.system(cloudy_run_command)

        sleep(1)

        if delete_out_file==True:
            os.remove(f'{run_name}.out')

        print('----------- Writing output files ------------- \n')


        with open(f'{run_name}_col_density.txt') as f:
            data=f.read()
            data=data.replace('#column density ','')

        with open(f'{run_name}_col_density.txt','w') as f:
            f.write(data)


        col_density=ascii.read(col_density)

        if n>=1:
            
            grid=loadtxt(grid_file,dtype=str)
            param_data=[]

            for i in range(n):
                param_data.append(grid[:,6+i].astype(float))

            for i in range(n):
                col_density.add_column(param_data[i],name=grid_parameters[i])
                
            param_data=tuple(zip(*param_data))

        if save_temp==True:

            data_temp=genfromtxt(temp_file,delimiter=[11,11,10,10,10])

            Te=data_temp[:,1]
            d2t_dr2=data_temp[:,4]

            log_Te=zeros(len(col_density))
            k=0

            for i,j in enumerate(d2t_dr2):
                if j==0:
                    log_Te[k]=round(log10(Te[i]),3)
                    k+=1

            col_density.add_column(log_Te,name='log_Te')

            if delete_temp_file==True:
                os.remove(temp_file)
                
        if n>=1:

            col_density.add_column(param_data,name='parameters')

        col_density.write(f'{run_name}_col_density_param.fits',overwrite=True)
        ascii.write(col_density,f'{run_name}_col_density_param.txt',format='ecsv',overwrite=True)
        col_density.write(f'../{run_name}_col_density_param.fits',overwrite=True)
            

        print('----------- Output files written -------------')
        os.chdir(cwd)

    else:

        os.system(f'rm -rf {path}')
        print('Cloudy run terminated...')


# run_name='component_III_PI_nH'

# qso='pg1116'
# z_abs=0.138527

# if not os.path.exists(f'{qso}/z={z_abs}'):
#     os.makedirs(f'{qso}/z={z_abs}')
    
# hden=[-5,1,0.02]
# metal=[-1]
# temp=None
# redshift=[0.138495,0.138506,0.138645]
# stop_nH=[14.97,13.6,16.04]

# ions=['H', 'H+', 'C+','C+2', 'C+3', 'N+', 'N+2', 'N+4', 'O','O+2','O+5','O+6','P+','Si+', 'Si+2', 'Si+3','Si+4','Fe+','Al+']

# for i in range(len(stop_nH)):

#     run_name=f'component_{toRoman(i+1)}_PI_nH'
#     run_cloudy(run_name, hden, metal, temp, redshift[i], stop_nH[i], ions, qso, z_abs, stop_T=50, save_temp=False, delete_temp_file=True, delete_out_file=True)

#     sleep(30)


'batch mode'


class ion():

    def __init__(self,name,z,b,logN):

        self.ion=name
        self.z=[x for x in z]
        self.b=[x for x in b]
        self.logN=[x for x in logN]
        self.comp=len(z)


class abs_system():

    def __init__(self,qso,z_abs,cont_mark='*'):

        file=f'../VPfit/{qso}/z={z_abs:.6f}/fit_params.txt'
        
        with open(file) as f:
            text=f.read()
            text=text.replace('A','')

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

        ions_all=unique(ions)
        ion_obj_dict={}

        for i in ions_all:
            mask=ions==i
            z_ion=z[mask]
            z_err_ion=z_err[mask]
            b_ion=b[mask]
            b_err_ion=b_err[mask]
            logN_ion=logN[mask]
            logN_err_ion=logN_err[mask]

            z_obj=[]
            b_obj=[]
            logN_obj=[]

            for j in range(len(z_ion)):
                z_obj.append([round(z_ion[j],6),round(z_err_ion[j],6)])
                b_obj.append([round(b_ion[j]),round(b_err_ion[j])])
                logN_obj.append([round(logN_ion[j],2),round(logN_err_ion[j],2)])

            obj=ion(i,z=z_obj,b=b_obj,logN=logN_obj)
            ion_obj_dict[i]=obj

        self.qso=qso
        self.z_abs=z_abs
        self.ion_obj=ion_obj_dict
        self.ions=ions_all[ions_all!='HI']
        self.n_ions=len(self.ions)


absorbers=[
            abs_system('he0056',0.043265),
            abs_system('pg1216',0.006328),
            abs_system('3c263',0.063397),
            abs_system('pg1222',0.054479),
            abs_system('rxj0439',0.005568),                        
            abs_system('uks0242',0.063850),
            abs_system('pg1259',0.046284),
            abs_system('pks1302',0.094839),
            abs_system('3c57',0.077430),
            abs_system('p1103',0.003934),
            abs_system('phl1811',0.080928),
           ]

n=0

hden=[-5,-3,1]
metal=[-1]
temp=None

ions=['H', 'H+', 'C+','C+2', 'C+3', 'N+', 'N+2', 'N+4', 'O','O+2','O+5','O+6','P+','Si+', 'Si+2', 'Si+3','Si+4','Fe+','Al+']


for a in absorbers[:2]:

    qso=a.qso
    z_abs=a.z_abs
    redshift=[x[0] for x in a.ion_obj['HI'].z]
    stop_nH=[x[0] for x in a.ion_obj['HI'].logN]

    if not os.path.exists(f'{qso}/z={z_abs}/logZ={metal[0]}'):
        os.makedirs(f'{qso}/z={z_abs}/logZ={metal[0]}')
        
  
    for i in range(len(stop_nH)):

        run_name=f'component_{toRoman(i+1)}_PI_nH'
        run_cloudy(run_name, hden, metal, temp, redshift[i], stop_nH[i], ions, qso, z_abs, stop_T=50, save_temp=True, delete_temp_file=True, delete_out_file=True)

        sleep(30)




































# H	H+	C+	C+0	C+0	N+0	N+0	O	O+0	O+0	O+0	Si+	Si+0	Si+0	Si+0
# 0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00
# 000000000	F	F	                  ok	0	000	-0.000000	-0.000000	-0.000000, -0.000000