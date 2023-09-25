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

for file in files:
    
    ind=[]
    for i,a in enumerate(file):
        if a=='_':
            ind.append(i)
        
    qso=file[ind[3]+1:ind[4]]

    data=ascii.read(f'Data/IGM_systems/{file}')

    line=data[col_name_dict['LINE_ID']]
    b=data[col_name_dict['BVALUE']]
    metal_lines=data[col_name_dict['NUM_METAL_LINES']]

    mask=logical_and(line=='Lya 1215',b>=40)
    mask=logical_and(mask,metal_lines>=3)

    data_BLA=data[mask]
    z_sys=data_BLA[col_name_dict['Z_SYS']]
    n=len(data_BLA)
    print(f'{qso} : {n} : {array(z_sys)}')
    k=k+n


print('\n ---------------------- \n')
print(k)