from astropy.table import Table,vstack
from numpy import *
from astropy.io import ascii,fits



def write_param(qso,z_abs,run_name,parameters,temp=False):

    # for i in range(0):

        # run_name = f'component_II_PI_nH_Z'
        # parameters=['log_nH','log_Z']
        print(run_name)

        with open(f'{qso}/z={z_abs}/{run_name}/{run_name}_col_density.txt') as f:
                data=f.read()
                data=data.replace('#column density ','')

        with open(f'{qso}/z={z_abs}/{run_name}/{run_name}_col_density.txt','w') as f:
            f.write(data)


        grid=loadtxt(f'{qso}/z={z_abs}/{run_name}/{run_name}_grid.txt',dtype=str)
        col_density=ascii.read(f'{qso}/z={z_abs}/{run_name}/{run_name}_col_density.txt')
        n=len(parameters)

        param_data=[]

        for i in range(n):
            param_data.append(grid[:,6+i].astype(float))

        for i in range(n):
            col_density.add_column(param_data[i],name=parameters[i])
            
        param_data=tuple(zip(*param_data))

        if temp==True:

            temp_file=f'{run_name}_temp.txt'
            
            data_temp=genfromtxt(f'{qso}/z={z_abs}/{run_name}/{temp_file}',delimiter=[11,11,10,10,10])

            Te=data_temp[:,1]
            d2t_dr2=data_temp[:,4]
            

            log_Te=zeros(len(col_density))
            # log_Te=[]
            # print(len(col_density))
            k=0
            
            for i,j in enumerate(d2t_dr2):
                
                if j==0:
                    print(k,Te[i],i)
                    log_Te[k]=round(log10(Te[i]),3)
                
                    k+=1
            
            # quit()  
            col_density.add_column(log_Te,name='log_Te')
            
        col_density.add_column(param_data,name='parameters')
        print(col_density.colnames,'\n')
        col_density.write(f'{qso}/z={z_abs}/{run_name}/{run_name}_col_density_param.fits',overwrite=True)
        ascii.write(col_density,f'{qso}/z={z_abs}/{run_name}/{run_name}_col_density_param.txt',format='ecsv',overwrite=True)
        col_density.write(f'{qso}/z={z_abs}/{run_name}_col_density_param.fits',overwrite=True)


def join_data():

    hdu=fits.open('Data/component_III_2d_sequence/component_III_2d_1_col_density_param.fits')
    data=Table(hdu[1].data)
    n=[len(data)]
    mask=data['H']!=0
    data=data[mask]
    n=[len(data)]


    for i in range(2,11):
        
        hdu1=fits.open(f'Data/component_III_2d_sequence/component_III_2d_{i}_col_density_param.fits')
        data1=Table(hdu1[1].data)
        n.append(len(data1))
        mask1=data1['H']!=0
        data1=data1[mask1]
        # n.append(len(data1))
        data=vstack([data,data1])

    print(n)
    print(sum(n))
    data.write(f'Data/component_III_2d_sequence/component_III_2d_joined_col_density_param.fits',overwrite=True)
    data.write(f'Data/component_III_2d_joined_col_density_param.fits',overwrite=True)

# join_data()

# qso='pks0637'
# z_abs=0.161064
# run_name = f'component_II_PI_nH'
# parameters=['log_nH']

# write_param(qso,z_abs,run_name,parameters,temp=True)


parameters=['log_nH','log_Z']
grid=loadtxt(f'pg0003/z=0.347579/component_II_CI/component_II_nH_Z_const_T_grid.txt',dtype=str)
col_density=ascii.read(f'pg0003/z=0.347579/component_II_CI/component_II_nH_Z_const_T_col_density.txt')
n=len(parameters)

param_data=[]

for i in range(n):
    param_data.append(grid[:,6+i].astype(float))

for i in range(n):
    col_density.add_column(param_data[i],name=parameters[i])
    
param_data=tuple(zip(*param_data))

# temp_file=f'scaling_approx/scaling_temp.txt'
            
# data_temp=genfromtxt(temp_file,delimiter=[11,11,10,10,10])

# Te=data_temp[:,1]
# d2t_dr2=data_temp[:,4]


# log_Te=zeros(len(col_density))
# k=0

# for i,j in enumerate(d2t_dr2):
    
#     if j==0:
#         # print(k,Te[i],i)
#         log_Te[k]=round(log10(Te[i]),3)
    
#         k+=1

# col_density.add_column(log_Te,name='log_Te')
             
col_density.add_column(param_data,name='parameters')
print(col_density.colnames,'\n')
col_density.write(f'pg0003/z=0.347579/component_II_CI/component_II_nH_Z_const_T_col_density_param.fits',overwrite=True)
ascii.write(col_density,f'pg0003/z=0.347579/component_II_CI/component_II_nH_Z_const_T_col_density_param.txt',format='ecsv',overwrite=True)

