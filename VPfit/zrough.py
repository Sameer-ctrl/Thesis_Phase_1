from astropy.io import ascii
from astropy.table import Table
import os
from numpy import argmin, linspace, logical_and, mean, sqrt, hstack, unique,vstack, where, zeros,array

spec='spec.fits'

data_system=ascii.read('Data/IGM_Danforth_Data/Systems/hlsp_igm_hst_cos_pg0003_g130m-g160m_v3_igm-systems.txt')
col_names_systems=['Z_SYS',  'DELTAV_SYS',  'WAVELENGTH',  'LINE_ID', 'z_ABS', 'SIGLEVEL', 'SNR', 'EQW', 'EQW_ERR', 'BVALUE', 'BVALUE_ERR', 'LOGN_lower_limit', 'LOGN', 'LOGN_ERR', 'FLAG_FIT', 'LCOD', 'NUM_SYS_LINES', 'NUM_METAL_LINES']
data_system.rename_columns(data_system.colnames,col_names_systems)

data_linelist=ascii.read('Data/IGM_Danforth_Data/Linelist/hlsp_igm_hst_cos_pg0003_g130m-g160m_v3_linelist.txt')
col_names_linelist=['WAVELENGTH','LINE_ID','z_ABS','SIGLEVEL','S/N','EQW','EQW_ERR','BVALUE','BVALUE_ERR','LOGN_lower_limit','LOGN','LOGN_ERR','FLAG_FIT','LCOD','FLAG_ISM','FLAG_AGN']
data_linelist.rename_columns(data_linelist.colnames,col_names_linelist)

z_sys=0.347
v_sep_lim=300

v_sys=3e5*(((1+z_sys)**2-1)/((1+z_sys)**2+1))

z1=sqrt((1+((v_sys-v_sep_lim)/3e5))/(1-((v_sys-v_sep_lim)/3e5)))-1
z2=sqrt((1+((v_sys+v_sep_lim)/3e5))/(1-((v_sys+v_sep_lim)/3e5)))-1

mask=logical_and(data_system['Z_SYS']>=z1, data_system['Z_SYS']<=z2)

data_system=data_system[mask]
z_sys=data_system['Z_SYS']
lines=data_system['LINE_ID'] 
z_lines=data_system['z_ABS']
b=data_system['BVALUE']
logN=data_system['LOGN']
wave=data_system['WAVELENGTH']


def init_guess_contamination(z_sys_uniq,ion):

    init_guess_list=[]
    table=[]

    for z in z_sys_uniq:

        mask=z_sys==z
      
        lines=data_system['LINE_ID'][mask]
        z_lines=data_system['z_ABS'][mask]
        b=data_system['BVALUE'][mask]
        logN=data_system['LOGN'][mask]
        wave=data_system['WAVELENGTH'][mask]
        ions=array([line.split(' ')[0] for line in lines])
        ions_rest_wave=array([int(line.split(' ')[1]) for line in lines])  

        mask=ions==ion
        guess_z=z_lines[mask]
        guess_b=b[mask]
        guess_logN=logN[mask]
        rest_wave=ions_rest_wave[mask]
        wave_range=[[w*(1+z1),w*(1+z2)] for w in rest_wave]

        init_guess=vstack((guess_z,zeros(len(guess_z)),guess_b,zeros(len(guess_z)),guess_logN,zeros(len(guess_z)))).transpose()

        if len(init_guess)>0:
            print(f'init_guess : {init_guess[0]}')
            print(f'range : {wave_range} \n')
            init_guess_list.append(init_guess[0])

    wave=data_linelist['WAVELENGTH']
    table=[]

    if len(init_guess)>0:

        for wr in wave_range:

            mask=logical_and(wave>=wr[0],wave<=wr[1])

            data=data_linelist[mask]
            line_contamination=data['LINE_ID']
            ions=array([line.split(' ')[0] for line in line_contamination])

            mask=ions!=ion
            data=data[mask]
            table.append(data)

        table=Table(vstack([*table]))  
    
    print(table)

    print(*init_guess_list)

uniq_z_sys=unique(z_sys)
init_guess_contamination(uniq_z_sys,'OVI')

quit()


class System():

    def __init__(self,z):

        mask=z_sys==z

        self.lines=data_system['LINE_ID'][mask]
        self.z_lines=data_system['z_ABS'][mask]
        self.b=data_system['BVALUE'][mask]
        self.logN=data_system['LOGN'][mask]
        self.wave=data_system['WAVELENGTH'][mask]
        self.ions=array([line.split(' ')[0] for line in self.lines])
        self.ions_rest_wave=array([int(line.split(' ')[1]) for line in self.lines])

    def param_guess(self,ion):

        mask=self.ions==ion
        guess_z=self.z_lines[mask]
        guess_b=self.b[mask]
        guess_logN=self.logN[mask]
        # wave=self.wave[mask]
        rest_wave=self.ions_rest_wave[mask]
        wave_range=[[w*(1+z1),w*(1+z2)] for w in rest_wave]

        init_guess=vstack((guess_z,zeros(len(guess_z)),guess_b,zeros(len(guess_z)),guess_logN,zeros(len(guess_z)))).transpose()

        if len(init_guess)>0:
            print(f'init_guess : {init_guess[0]}')
            print(f'range : {wave_range} \n')


            wave=data_linelist['WAVELENGTH']
            table=[]

            for wr in wave_range:

                mask=logical_and(wave>=wr[0],wave<=wr[1])

                data=data_linelist[mask]
                line_contamination=data['LINE_ID']
                ions=array([line.split(' ')[0] for line in line_contamination])

                mask=ions!=ion
                data=data[mask]
                table.append(data)

            table=Table(vstack([*table]))  
            
            print(table)
            print('\n -------------------------------------------------------------------- \n')   

uniq_z_sys=unique(z_sys)

# sys1=System(uniq_z_sys[1])  
# sys1.param_guess('OVI')  
# sys1=System(uniq_z_sys[2])  
# sys1.param_guess('OVI') 

for z in uniq_z_sys:
    sys=System(z)  
    sys.param_guess('SiII') 


