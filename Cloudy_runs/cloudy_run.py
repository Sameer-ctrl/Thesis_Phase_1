import os
from numpy import *
from astropy.io import ascii
from time import sleep


def run_cloudy(run_name, hden, metal, temp, redshift, stop_nH, ions, stop_T=1000, save_temp=False, save_Hyd=False, uvb_scale = 1, uvb_Q = 18, He_abun=0.08156498, delete_out_file=False, delete_temp_file=True):

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

    lines=[uv_b,abundance,hden_line,metal_line,temp_line,stop_criteria_nH,save_grid,save_Temp,save_hyd,save_col_den]

    os.mkdir(f'Data/{run_name}')
    file=open(f'Data/{run_name}/{run_name}.in','w+')
    file.writelines(lines)
    file.write('end')
    file.close()


    file=open(f'Data/{run_name}/{run_name}.in','r')

    print('-------------Cloudy input file------------- \n')
    print(file.read())
    print('\n--------------------------------------- \n')

    if n>=1:
        print(f'Grid parameters : {grid_parameters}')
        # print(f'No. of models : {n_models} \n')

    a=input('Check the input file for cloudy. \nPress Y to continue.\n')

    if a=='y':

        os.chdir(f'Data/{run_name}')

        cloudy_path_desktop='/home/sameer/Documents/Sameer/c22.02/source/cloudy.exe'
        cloudy_path_workstation='/home/sameer/Sameer/c22.02/source/cloudy.exe'
        cloudy_path_pc='/home/sameer/cloudy/source/sys_gcc/cloudy.exe'

        cloudy_path=cloudy_path_desktop

        cloudy_run_command=f'{cloudy_path} -r {run_name}'

        print(f'\n Cloudy running in Data/{run_name} ... \n')
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

            if delete_temp_file==True:
                os.remove(temp_file)
                
        if n>=1:

            col_density.add_column(param_data,name='parameters')

        col_density.write(f'{run_name}_col_density_param.fits',overwrite=True)
        ascii.write(col_density,f'{run_name}_col_density_param.txt',format='ecsv',overwrite=True)
        col_density.write(f'../{run_name}_col_density_param.fits',overwrite=True)
            

        print('----------- Output files written -------------')

    else:

        # os.remove(f'/home/sameer/Sameer/Thesis_Phase_1/Cloudy_runs/{run_name}')
        print('Cloudy run terminated...')


run_name='component_II_nH_const_T'

hden=[-5,0,1]
metal=[-1]
temp=[5.29]
redshift=0.34758
stop_nH=14.13

ions=['H', 'H+', 'C+','C+2', 'C+3', 'O+5','O+6','Si+', 'Si+2', 'Si+3']


run_cloudy(run_name, hden, metal, temp, redshift, stop_nH, ions)

