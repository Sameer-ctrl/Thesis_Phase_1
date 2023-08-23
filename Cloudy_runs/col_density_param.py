from astropy.table import Table,vstack
from numpy import loadtxt,sum, arange
from astropy.io import ascii,fits


def write_param():

    # for i in range(0):

        run_name = f'component_III_nH_Z'
        parameters=['log_nH','log_Z']
        print(run_name)

        # with open(f'{run_name}/{run_name}_col_density.txt') as f:
        #         data=f.read()
        #         data=data.replace('#column density ','')

        # with open(f'{run_name}/{run_name}_col_density.txt','w') as f:
        #     f.write(data)


        grid=loadtxt(f'{run_name}/{run_name}_grid.txt',dtype=str)
        col_density=ascii.read(f'{run_name}/{run_name}_col_density.txt')
        n=len(parameters)

        param_data=[]

        for i in range(n):
            param_data.append(grid[:,6+i].astype(float))

        for i in range(n):
            col_density.add_column(param_data[i],name=parameters[i])
            
        param_data=tuple(zip(*param_data))

        col_density.add_column(param_data,name='parameters')
        print(col_density.colnames,'\n')
        col_density.write(f'{run_name}/{run_name}_col_density_param.fits',overwrite=True)
        ascii.write(col_density,f'{run_name}/{run_name}_col_density_param.txt',format='ecsv',overwrite=True)
        col_density.write(f'Data/{run_name}_col_density_param.fits',overwrite=True)


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
write_param()