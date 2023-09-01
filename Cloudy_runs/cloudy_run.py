import os
from numpy import *
from astropy.io import ascii
from time import sleep

run_name = 'component_II_T1'
parameters=['log_T']

'UV background model : KS19'

redshift = 0.34758
uvb_scale = 1
uvb_Q = 18

He_abun=0.08156498

hden = [-4]
metal = [0] 
# print((len(arange(*hden))+1)*(len(arange(*metal))+1))
# quit()
# print((len(arange(*hden))+1))
temp = [4,7,0.02]
# print((len(arange(*hden))+1)*(len(arange(*temp))+1))


stop_nH=14.13
stop_T=20

grid_file=f'{run_name}_grid.txt'
# hyd_cond=f'{run_name}_hydrogen.txt'
col_density=f'{run_name}_col_density.txt'
temp_file=f'{run_name}_temp.txt'

ions=['H', 'H+', 'C+','C+2', 'C+3', 'O+5','O+6','Si+', 'Si+2', 'Si+3']

uv_b=f'TABLE KS19 redshift = {redshift} scale = {uvb_scale} Q = {uvb_Q} \n'
abundance=f'abundances "solar_GASS10.abn" \nelement helium abundance {He_abun} linear \n'
stop_criteria_nH=f'stop column density {stop_nH} neutral H \nstop temperature {stop_T} K linear \n'
save_grid=f'save grid "{grid_file}" last no hash \n'
# save_hyd=f'save hydrogen conditions "{hyd_cond}" last no clobber \n'
save_col_den=f'save species column density "{col_density}" no hash \n'
save_temp=f'save temperature "{temp_file}" no hash \n'

for i in ions:
    save_col_den+=f'"{i}" \n'

if len(hden)>1:
    hden = f'hden -5 log vary \ngrid range from {hden[0]} to {hden[1]} with {hden[2]} dex step \n'

else:
    hden = f'hden {hden[0]} log \n'

if len(metal)>1:
    metal = f'metals -1 log vary \ngrid range from {metal[0]} to {metal[1]} with {metal[2]} dex step \n'

else:
    metal = f'metals {metal[0]} log \n'

if temp!=None:

    if len(temp)>1:
        temp = f'constant temperature 5 log vary \ngrid range from {temp[0]} to {temp[1]} with {temp[2]} dex step \n'

    else:
        temp = f'constant temperature {temp[0]} log \n'

else:
    temp=''

lines=[uv_b,abundance,hden,metal,temp,stop_criteria_nH,save_grid,save_col_den]

os.mkdir(run_name)
file=open(f'{run_name}/{run_name}.in','w+')
file.writelines(lines)
file.write('end')
file.close()


file=open(f'{run_name}/{run_name}.in','r')

print('-------------Cloudy input file------------- \n')
print(file.read())
print('\n--------------------------------------- \n')
print('\n--------- Have you checked varying parameters ? ----------')
print(parameters,'\n')

a=input('Check the input file for cloudy. \n Press Y to continue.')

if a=='y':

    os.chdir(run_name)

    cloudy_path='/home/sameer/Sameer/c22.02/source/cloudy.exe'
    cloudy_run_command=f'{cloudy_path} -r {run_name}'

    print(f'\n Cloudy running in {run_name} ... \n')
    os.system(cloudy_run_command)

    sleep(1)


    print('----------- Writing output files ------------- \n')


    with open(f'{run_name}_col_density.txt') as f:
        data=f.read()
        data=data.replace('#column density ','')

    with open(f'{run_name}_col_density.txt','w') as f:
        f.write(data)


    grid=loadtxt(grid_file,dtype=str)
    col_density=ascii.read(col_density)

    n=len(parameters)

    param_data=[]

    for i in range(n):
        param_data.append(grid[:,6+i].astype(float))

    for i in range(n):
        col_density.add_column(param_data[i],name=parameters[i])
        
    param_data=tuple(zip(*param_data))

    # data_temp=genfromtxt(temp_file,delimiter=[11,11,9,10,10])

    # Te=data_temp[:,1]
    # d2t_dr2=data_temp[:,4]

    # log_Te=zeros(len(col_density))
    # k=0

    # for i,j in enumerate(d2t_dr2):
    #     if j==0:
    #         log_Te[k]=round(log10(Te[i]),3)
    #         k+=1

    # col_density.add_column(log_Te,name='log_Te')
    col_density.add_column(param_data,name='parameters')


    col_density.write(f'{run_name}_col_density_param.fits',overwrite=True)
    ascii.write(col_density,f'{run_name}_col_density_param.txt',format='ecsv',overwrite=True)
    col_density.write(f'../Data/{run_name}_col_density_param.fits',overwrite=True)

    print('----------- Output files written -------------')

else:

    # os.remove(f'/home/sameer/Sameer/Thesis_Phase_1/Cloudy_runs/{run_name}')
    print('Cloudy run terminated...')
    
