from astropy.io import ascii
import os
from numpy import *
import matplotlib.pyplot as plt
import matplotlib

# matplotlib.rcParams['backend']='TkAgg'
# quit()

v_z=lambda z : 3e5*(((1+round(z,6))**2-1)/((1+round(z,6))**2+1))  # v at z
z_v=lambda v : sqrt((1+((v)/3e5))/(1-((v)/3e5)))-1      # z at v


col_names=['Z_SYS',  'DELTAV_SYS',  'WAVELENGTH',  'LINE_ID', 'z_ABS', 'SIGLEVEL', 'SNR', 'EQW', 'EQW_ERR', 'BVALUE', 'BVALUE_ERR', 'LOGN_lower_limit', 'LOGN', 'LOGN_ERR', 'FLAG_FIT', 'LCOD', 'NUM_SYS_LINES', 'NUM_METAL_LINES']

n_BLA_metal_lines=0        # systems with BLA and >=3 metal lines
n_BLA_metal_ions=0         # systems with BLA and >=3 distinct metal ions
n_los_BLA_ions=0           # LOS with BLA and >=3 distinct metal ions
n_los_BLA_OVI_ions=0       # LOS with BLA, OVI and >=2 other distinct metal ions
n_los_BLA_OVI=0            # LOS with BLA and OVI
n_los_BLA=0                # LOS with BLA
n_BLA=0                    # systems with BLA
n_BLA_OVI=0                # systems with BLA and OVI
n_BLA_OVI_metal_ions=0     # systems with BLA, OVI and >=2 other distinct metal ions

'''

Class I     : BLA and OVI
Class II    : BLA and ditinct metal ions >=3
Class III   : BLA, OVI and other metal ions >=2

'''

files=os.listdir('Data/IGM_Danforth_Data/Systems')
# qso_OVI=unique(['3c263', 'pks0637', 'pks0637', 'pg1424', 'pg0003', 'pg0003', 'pg0003', 'pg1216', 's135712', '1es1553', 'sbs1108', 'pg1222', 'pg1116', 'h1821', 'h1821', 'pg1121', 'pks0405'])


for file in files:
# for q in qso_OVI:
        
    qso=file[:-16]
    print(f'\n{qso}')

    data=ascii.read(f'Data/IGM_Danforth_Data/Systems/{file}')
    data.rename_columns(data.colnames,col_names)

    lines=data['LINE_ID']
    z_sys=data['Z_SYS']
    b=data['BVALUE']
    n_metal_lines=data['NUM_METAL_LINES']

    'masking systems by BLA'

    mask_BLA=logical_and(lines=='Lya 1215',b>=45)
    z_sys_BLA=data['Z_SYS'][mask_BLA]
    n_BLA+=len(z_sys_BLA)

    if len(z_sys_BLA)>0:
        n_los_BLA+=1

    z_plot_BLA_OVI=[]
    a=[]  # to count BLA-OVI los
    b=[]  # to count BLA-OVI and ions>2 los

    for z in z_sys_BLA:
            
        'masking systems having BLA by z'


        v_z1=v_z(z)
        z1=z_v(v_z1-500)
        z2=z_v(v_z1+500)

        # mask_sys=z_sys==z
        mask_sys=logical_and(z_sys>z1, z_sys<z2)
        line_sys=lines[mask_sys].value
        metal_lines=[]

        for line in line_sys:
            if line[:2]!='Ly':
                metal_lines.append(line)
        
        ions=set([m_line.split(' ')[0] for m_line in metal_lines])

        if 'OVI' in ions:

            a.append(True)
            
            n_BLA_OVI+=1
            z_plot_BLA_OVI.append(z)
            # print(f'BLA and OVI  : {z:.6f} : {metal_lines} : {ions} : {len(ions)}')
        


    mask_BLA_metal_lines=logical_and(mask_BLA,n_metal_lines>=3)
    z_sys_BLA_metal_lines=data['Z_SYS'][mask_BLA_metal_lines]
    n=len(z_sys_BLA_metal_lines)
    n_BLA_metal_lines+=n

    if n>0:

        i=0
        z_plot_BLA_OVI_metal_ions=[]
        for z in z_sys_BLA_metal_lines:
            
            'maksing systems having BLA and metal lines >=3 by z'

            v_z1=v_z(z)
            z1=z_v(v_z1-500)
            z2=z_v(v_z1+500)

            # mask_sys=z_sys==z
            mask_sys=logical_and(z_sys>z1, z_sys<z2)

            line_sys=lines[mask_sys].value
            metal_lines=[]

            for line in line_sys:
                if line[:2]!='Ly':
                    metal_lines.append(line)
            
            ions=set([m_line.split(' ')[0] for m_line in metal_lines])
            
            'systems with distinct metal ions >= 3'

            if len(ions)>=3:
                n_BLA_metal_ions+=1
                i=1

                # print(f'Ions > 2     : {z:.6f} : {metal_lines} : {ions} : {len(ions)}')

                'systems with OVI and other distinct metal ions >= 2'

                if 'OVI' in ions:
                    n_BLA_OVI_metal_ions+=1
                    b.append(True)
                    z_plot_BLA_OVI_metal_ions.append(z)
                    print(f'OVI and ions : {z:.6f} : {metal_lines} : {ions} : {len(ions)}')

        n_los_BLA_ions+=i

        # plt.title(f'LOS : {qso}')
        # plt.scatter(z_plot_BLA_OVI,ones(len(z_plot_BLA_OVI)),label='BLA and OVI')
        # plt.scatter(z_plot_BLA_OVI_metal_ions,2*ones(len(z_plot_BLA_OVI_metal_ions)),label='BLA, OVI and metal ions')
        # plt.legend()
        # plt.show()

    if any(a):
        n_los_BLA_OVI+=1
    
    if any(b):
        n_los_BLA_OVI_ions+=1

        

print('\n----------------------------------\n')

print(f'Systems with BLA and metal lines >=3               : {n_BLA_metal_lines}')
print(f'Systems with BLA and distinct metal ions >=3       : {n_BLA_metal_ions}')
print(f'No. of los with BLA                                : {n_los_BLA}')
print(f'No. of los with BLA and OVI                        : {n_los_BLA_OVI}')
print(f'No. of los BLA, OVI and other metal ions >=2       : {n_los_BLA_OVI_ions}')
print(f'No. of los with BLA and distinct metal ions >= 3   : {n_los_BLA_ions}')
print(f'Systems with BLA                                   : {n_BLA}')
print(f'Systems with BLA and OVI                           : {n_BLA_OVI}')
print(f'Systems with BLA, OVI and other metal ions >=2     : {n_BLA_OVI_metal_ions}')

