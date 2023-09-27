from astropy.io import ascii
from astropy.table import Table
import os
from numpy import argmin, linspace, logical_and, mean, sqrt, hstack, unique,vstack, where, zeros,array


spec='spec.fits'

data_system=ascii.read('Data/IGM_Danforth_Data/Systems/hlsp_igm_hst_cos_pg1121_g130m-g160m_v3_igm-systems.txt')
col_names_systems=['Z_SYS',  'DELTAV_SYS',  'WAVELENGTH',  'LINE_ID', 'z_ABS', 'SIGLEVEL', 'SNR', 'EQW', 'EQW_ERR', 'BVALUE', 'BVALUE_ERR', 'LOGN_lower_limit', 'LOGN', 'LOGN_ERR', 'FLAG_FIT', 'LCOD', 'NUM_SYS_LINES', 'NUM_METAL_LINES']
data_system.rename_columns(data_system.colnames,col_names_systems)

data_linelist=ascii.read('Data/IGM_Danforth_Data/Linelist/hlsp_igm_hst_cos_pg1121_g130m-g160m_v3_linelist.txt')
col_names_linelist=['WAVELENGTH','LINE_ID','z_ABS','SIGLEVEL','S/N','EQW','EQW_ERR','BVALUE','BVALUE_ERR','LOGN_lower_limit','LOGN','LOGN_ERR','FLAG_FIT','LCOD','FLAG_ISM','FLAG_AGN']
data_linelist.rename_columns(data_linelist.colnames,col_names_linelist)

lsf_files=os.listdir('Data/COS_LSF')

def lsf_file_search(lambda1,lambda2):

    cen_wave=mean([lambda1,lambda2])
    wave_diff=[abs(cen_wave-int(file[14:18])) for file in lsf_files]
    i=argmin(wave_diff)

    return lsf_files[i]

z_sys=0.192
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



quit()  











ion_to_fit='OVI'

for z in z_sys:
    mask=z_sys==z
    line_sys=lines[mask]

    ions=[line.split(' ')[0] for line in line_sys]

    ind=[]
    for i,ion in enumerate(ions):
        if ion==ion_to_fit:
            ind.append(i)

    guess_z=z_lines[ind[0]]
    guess_b=b[ind[0]]
    guess_logN=logN[ind[0]]

    init_guess=vstack((guess_z,0,guess_b,0,guess_logN,0)).transpose()
    # init_guess=vstack((guess_z,zeros(len(ind)),guess_b,zeros(len(ind)),guess_logN,zeros(len(ind)))).transpose()
    print(z)
    print(init_guess)
    print('------------------------')






