from numpy import *
from astropy.io import ascii
import matplotlib.pyplot as plt
from scipy.constants import c,e,m_e


data=ascii.read('Data/PG0003+158_rebinned_cont_norm.asc')
atom_data=loadtxt('Data/atom.txt',dtype=str,usecols=(0,1,2))

ion=atom_data[:,0]
rest_wave=atom_data[:,1].astype(float)
f=atom_data[:,2]

line=[]

for i in range(len(atom_data)):
    line.append(f'{ion[i]}_{int(rest_wave[i])}')

wave=data['WAVE']
flux=data['FLUX']

aod=-log(flux)

ions={'OVI_1031':[1391.6,1391.4],'OVI_1037':[1398.3,1398.9],'CII_1036':[1396.58,1397],'CIII_977':[1316.5 1317.4],'SiII_1260':[1698.75,1699.05],'SiIII_1206':[1626.05,1626.5]}
# wave_range=[]
'''
SiII: 1698.75,1699.05
SiIII: 1626.05,1626.5
CII: 1396.58,1397
CIII: 1316.5 1317.4
OVI_1031 : 1391.6,1391.4
OVI_1037 : 1398.3,1398.9

'''


def osc_strength(ions_dict):

    ions=ions_dict.keys()
    osc=[]

    for ion in ions:
        j=0
        for i in line:

            if i==ion:
                break
            else:
                j+=1

        osc.append(float(f[j]))
    
    return osc








