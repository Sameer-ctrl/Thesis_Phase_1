from astropy.io import ascii
import os
from numpy import *


col_names=['Z_SYS',  'DELTAV_SYS',  'WAVELENGTH',  'LINE_ID', 'z_ABS', 'SIGLEVEL', 'SNR', 'EQW', 'EQW_ERR', 'BVALUE', 'BVALUE_ERR', 'LOGN_lower_limit', 'LOGN', 'LOGN_ERR', 'FLAG_FIT', 'LCOD', 'NUM_SYS_LINES', 'NUM_METAL_LINES']
col_num=['col1', 'col2', 'col3', 'col4', 'col5', 'col6', 'col7', 'col8', 'col9', 'col10', 'col11', 'col12', 'col13', 'col14', 'col15', 'col16', 'col17', 'col18']

col_name_dict=dict(zip(col_names,col_num))

files=os.listdir('Data/IGM_systems')


for file in files:

    ind=[]
    for i,a in enumerate(file):
        if a=='_':
            ind.append(i)
        
    qso=file[ind[3]+1:ind[4]]

k=0
j=0
n_qso=0

for file in files:
    
    ind=[]
    for i,a in enumerate(file):
        if a=='_':
            ind.append(i)
        
    qso=file[ind[3]+1:ind[4]]

    data=ascii.read(f'Data/IGM_systems/{file}')

    line=data[col_name_dict['LINE_ID']]
    z_sys=data[col_name_dict['Z_SYS']]
    b=data[col_name_dict['BVALUE']]
    n_metal_lines=data[col_name_dict['NUM_METAL_LINES']]

    mask=logical_and(line=='Lya 1215',b>=40)
    mask=logical_and(mask,n_metal_lines>=3)

    data_BLA=data[mask]
    n=len(data_BLA)

    if n>0:

        z_sys_BLA=data_BLA[col_name_dict['Z_SYS']]
        print(f'\n {qso} : {n}')
        m=0

        for z in z_sys_BLA:

            mask=z_sys==z
            line_sys=line[mask].value
            metal_lines=[]

            for l in line_sys:
                if l[:2]!='Ly':
                    metal_lines.append(l)
            
            ions=set([m.split(' ')[0] for m in metal_lines])
            if len(ions)>=3:
                m=m+1

            print(f'{z} : {metal_lines} : {ions} :{len(ions)}')

        j=j+m
        n_qso=n_qso+1

    k=k+n
        

print('\n ---------------------- \n')
print(k)
print(j)
print(n_qso)