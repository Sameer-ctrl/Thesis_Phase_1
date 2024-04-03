from operator import le
from astropy.io import ascii
from astropy.table import Table
import os
import subprocess
from numpy import *
import matplotlib.pyplot as plt


lsf_files=os.listdir('../VPfit/Data/COS_LSF')

qso_list=['pg1424', 'sbs1108', 'pks0405', 'pks0637', 'pks0637', 'pks0637', 'pg0003', 'pg0003', 's135712', 'pg1216', '1es1553', 'pg1222', 'h1821', 'h1821', '3c263', 'pg1121', 'pg1116']
z_abs_list=[0.146789, 0.463201, 0.167125, 0.161068, 0.161068, 0.417573, 0.386094, 0.42188, 0.097767, 0.282195, 0.187731, 0.3786, 0.170062, 0.224832, 0.140754, 0.192434, 0.138527]


for i in range(len(qso_list)):

    qso=qso_list[i]
    z_abs=z_abs_list[i]

    print(f'{qso} : {z_abs}')

    data_system=ascii.read(f'../VPfit/Data/IGM_Danforth_Data/Systems/{qso}_igm-systems.txt')
    col_names_systems=['Z_SYS',  'DELTAV_SYS',  'WAVELENGTH',  'LINE_ID', 'z_ABS', 'SIGLEVEL', 'SNR', 'EQW', 'EQW_ERR', 'BVALUE', 'BVALUE_ERR', 'LOGN_lower_limit', 'LOGN', 'LOGN_ERR', 'FLAG_FIT', 'LCOD', 'NUM_SYS_LINES', 'NUM_METAL_LINES']
    data_system.rename_columns(data_system.colnames,col_names_systems)

    data_linelist=ascii.read(f'../VPfit/Data/IGM_Danforth_Data/Linelist/{qso}_linelist.txt')
    col_names_linelist=['WAVELENGTH','LINE_ID','z_ABS','SIGLEVEL','S/N','EQW','EQW_ERR','BVALUE','BVALUE_ERR','LOGN_lower_limit','LOGN','LOGN_ERR','FLAG_FIT','LCOD','FLAG_ISM','FLAG_AGN']
    data_linelist.rename_columns(data_linelist.colnames,col_names_linelist)

    z_sys=data_system['Z_SYS']
    lines_sys=data_system['LINE_ID']

    mask1=z_sys==z_abs
    mask2=logical_or(lines_sys=='Lya 1215',lines_sys=='OVI 1032')

    mask=logical_and(mask1,mask2)

    data_system=data_system[mask]
    z_sys=data_system['Z_SYS']
    lines=data_system['LINE_ID'] 

    print(data_system.pprint_all())
    print('\n')

quit()

ion_dict={'HI':'Lya 1215','OVI':'OVI 1032'}

def param_print(qso,z_abs):

    print(f'{qso} : {z_abs}')

    for ion in ['HI','OVI']:

        # spec=f'{qso}_cont_norm.asc'
        v_sep_lim=300


        data_system=ascii.read(f'../VPfit/Data/IGM_Danforth_Data/Systems/{qso}_igm-systems.txt')
        col_names_systems=['Z_SYS',  'DELTAV_SYS',  'WAVELENGTH',  'LINE_ID', 'z_ABS', 'SIGLEVEL', 'SNR', 'EQW', 'EQW_ERR', 'BVALUE', 'BVALUE_ERR', 'LOGN_lower_limit', 'LOGN', 'LOGN_ERR', 'FLAG_FIT', 'LCOD', 'NUM_SYS_LINES', 'NUM_METAL_LINES']
        data_system.rename_columns(data_system.colnames,col_names_systems)

        data_linelist=ascii.read(f'../VPfit/Data/IGM_Danforth_Data/Linelist/{qso}_linelist.txt')
        col_names_linelist=['WAVELENGTH','LINE_ID','z_ABS','SIGLEVEL','S/N','EQW','EQW_ERR','BVALUE','BVALUE_ERR','LOGN_lower_limit','LOGN','LOGN_ERR','FLAG_FIT','LCOD','FLAG_ISM','FLAG_AGN']
        data_linelist.rename_columns(data_linelist.colnames,col_names_linelist)


        v_absorber=3e5*(((1+z_abs)**2-1)/((1+z_abs)**2+1))

        z1=sqrt((1+((v_absorber-v_sep_lim)/3e5))/(1-((v_absorber-v_sep_lim)/3e5)))-1
        z2=sqrt((1+((v_absorber+v_sep_lim)/3e5))/(1-((v_absorber+v_sep_lim)/3e5)))-1

        mask=logical_and(data_system['Z_SYS']>=z1, data_system['Z_SYS']<=z2)

        data_system=data_system[mask]
        z_sys=data_system['Z_SYS']
        lines=data_system['LINE_ID'] 
        ions_sys=array([line.split(' ')[0] for line in lines])

        def init_guess_contamination(z_sys_uniq,ion):

            init_guess_list=[]
            wave_range_list=[]
            table=[]

            for z in z_sys_uniq:

                if ion=='HI':

                    ions_sys_HI=array([i[:2] for i in ions_sys])
                    mask=logical_and(z_sys==z,ions_sys_HI=='Ly')

                else:
                    mask=logical_and(z_sys==z,ions_sys==ion)

                lines=data_system['LINE_ID'][mask]
                guess_z=data_system['z_ABS'][mask]
                guess_b=data_system['BVALUE'][mask]
                guess_logN=data_system['LOGN'][mask]
                wave=data_system['WAVELENGTH'][mask]
                ions_rest_wave=array([int(line.split(' ')[1]) for line in lines]) 
            
                wave_range=[[w*(1+z1),w*(1+z2)] for w in ions_rest_wave]

                init_guess=vstack((guess_z,zeros(len(guess_z)),guess_b,zeros(len(guess_z)),guess_logN,zeros(len(guess_z)))).transpose()

                if len(init_guess)>0:
                    # print(f'init_guess : {init_guess[0]}')
                    # print(f'range : {wave_range} \n')
                    init_guess_list.append(init_guess[0])

                    for wr in wave_range:
                        wave_range_list.append(wr)

            wave=data_linelist['WAVELENGTH']
            table=[]

            if len(init_guess_list)>0:

                for wr in wave_range_list:

                    mask=logical_and(wave>=wr[0],wave<=wr[1])

                    data=data_linelist[mask]
                    table.append(data)

                table=Table(hstack([*table]))  
                
            unique_wave_range=[]

            for wr in wave_range_list:
                if wr not in unique_wave_range:
                    unique_wave_range.append(wr)
            
            wave_range=vstack(unique_wave_range)
            init_guess=vstack(init_guess_list)

            # print('\n--------------------- See contamination if there any ------------------------\n')
            
            mask1=table['LINE_ID']==ion_dict[ion]

            if ion=='HI':
                mask2=table['BVALUE']>=45
                mask=logical_and(mask1,mask2)
                table[mask].pprint_all()

            else:
                table[mask1].pprint_all()

            print(f'\n----------------- {ion} ----------------------\n')
            print(init_guess)
            print('\n')
            # print('\n----------------- Wavelength Ranges --------------------\n')
            # print(wave_range)
            # print('----------------------------------------------------------\n')

            return wave_range,init_guess




        uniq_z_sys=unique(z_sys)
        wave_range,init_guess=init_guess_contamination(uniq_z_sys,ion)


