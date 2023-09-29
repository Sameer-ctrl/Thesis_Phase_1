from astropy.io import ascii
from astropy.table import Table
import os
from numpy import *


qso='pks0405'
spec=f'{qso}_cont_norm.asc'
z_absorber=0.167125
v_sep_lim=300
ion='CII'

lsf_files=os.listdir('Data/COS_LSF')

data_system=ascii.read(f'Data/IGM_Danforth_Data/Systems/{qso}_igm-systems.txt')
col_names_systems=['Z_SYS',  'DELTAV_SYS',  'WAVELENGTH',  'LINE_ID', 'z_ABS', 'SIGLEVEL', 'SNR', 'EQW', 'EQW_ERR', 'BVALUE', 'BVALUE_ERR', 'LOGN_lower_limit', 'LOGN', 'LOGN_ERR', 'FLAG_FIT', 'LCOD', 'NUM_SYS_LINES', 'NUM_METAL_LINES']
data_system.rename_columns(data_system.colnames,col_names_systems)

data_linelist=ascii.read(f'Data/IGM_Danforth_Data/Linelist/{qso}_linelist.txt')
col_names_linelist=['WAVELENGTH','LINE_ID','z_ABS','SIGLEVEL','S/N','EQW','EQW_ERR','BVALUE','BVALUE_ERR','LOGN_lower_limit','LOGN','LOGN_ERR','FLAG_FIT','LCOD','FLAG_ISM','FLAG_AGN']
data_linelist.rename_columns(data_linelist.colnames,col_names_linelist)


v_absorber=3e5*(((1+z_absorber)**2-1)/((1+z_absorber)**2+1))

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

    if len(init_guess)>0:

        for wr in wave_range:

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

    print('\n--------------------- See contamination if there any ------------------------\n')
    table.pprint_all()
    print('\n----------------- Initial guesses ----------------------\n')
    print(init_guess)
    print('\n----------------- Wavelength Ranges --------------------\n')
    print(wave_range)
    print('----------------------------------------------------------\n')

    return wave_range,init_guess


def lsf_file_search(lambda1,lambda2):

    cen_wave=mean([lambda1,lambda2])
    wave_diff=[abs(cen_wave-int(file[14:18])) for file in lsf_files]
    i=argmin(wave_diff)

    return lsf_files[i]

uniq_z_sys=unique(z_sys)
wave_range,init_guess=init_guess_contamination(uniq_z_sys,ion)

spectrum_lines=''

for i,wr in enumerate(wave_range):
    lsf_file=lsf_file_search(*wr)
    line=f'%% {spec} {i+1} {wr[0]:.3f} {wr[1]:.3f} pfinst=LSF/{lsf_file} !\n'
    spectrum_lines+=line

guess_line=''

for i in init_guess:

    i_str=''

    for param in i:
        i_str+=f'{param}   '

    line=f'{ion} {i_str} 0 ! \n'
    guess_line+=line


with open(f'{qso}/fit_{ion}.asc','w') as f:
    f.writelines([spectrum_lines,guess_line])

with open(f'{qso}/fit_{ion}.asc','r') as f:
    print(f.read())

