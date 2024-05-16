import os
from numpy import *
from astropy.io import ascii
import matplotlib.pyplot as plt
from scipy.integrate import simpson



z_abs=0.347586
file_path=f'pg0003/z={z_abs:.6f}/VPfit_chunks'

data=loadtxt('../Python/Data/rest_wave.txt',dtype=str)

line_atom=data[:,0]
wave_rest=data[:,1].astype(float)

rest_wave={}

for i in range(len(line_atom)):
    rest_wave.update({line_atom[i]:wave_rest[i]})


def cos_siglevel(W, wavelength, b, snpix=None, snx=None, disp=None, binning=None, xopt=None, fcx=None):

    if snpix is None and snx is None:
        print('SIGLEVEL: either SNPIX or SNX must be specified!!!')
        return -1

    if disp is None:

        if wavelength <= 1425:  #1425
            disp=9.97

        else:
            disp=12.23
    
    # if disp is None:

    #     if wavelength <= 1460:
    #         disp=29.9 

    #     else:
    #         disp=36.7

    if binning is None:
        binning = 1

    dx=(1e3*b*wavelength)/(2.998e5*disp)

    xopt=(1.605*dx)+(5.1*(dx**(-0.25)))
    eta=0.15+(xopt**0.37)
    fcx=0.743-(0.185*exp(-dx/11.6))


    if snx is not None:
        siglevel=snx*(W/disp)*(fcx/xopt)
    
    else:
        if binning==1:
            sn1=snpix 

        else:
            sn1=snpix/(0.15+(binning**0.37))

        siglevel=sn1*(W/disp)*eta*(fcx/xopt)

    # xopt/=binning

    return siglevel


def eq_w(line,line_name):

    data=loadtxt(f'{file_path}/{line}.txt',comments='!')

    cen_wave_rest=rest_wave[line_name]

    wave=data[:,0]
    cont=data[:,3]

    w=wave[-1]-wave[0]-simpson(cont,wave)
    
    return w*1000  #mA,


file=f'pg0003/z={z_abs:.6f}/fit_params.txt'
        
with open(file) as f:
    text=f.read()
    text=text.replace('A','')
    # print(text)
    # print('\n')

with open('temp_param_file.txt','w') as f:
    f.write(text)

param_file=loadtxt('temp_param_file.txt',dtype=str,comments=('>','#'))

ions=param_file[:,0]
mask=[]

for i in ions:
    mask.append('*' not in i)  

ions=ions[mask]
z=param_file[:,1].astype(float)[mask]
b=param_file[:,3].astype(float)[mask]

cen_wave=(1+z)*array([1260.4221,1206.5,1036.3367,977.0201,1031.927,1031.927,1215.6701,1215.6701,1215.6701])

snpix=[12,11,15,16,18,19,12,12,12]

lines=['SiII_1260', 'SiIII_1206', 'CII_1036', 'CIII_977_1', 'OVI_1032_1', 'OVI_1032_2', 'HI_1215_1', 'HI_1215_2', 'HI_1215_3']
line_name=['SiII_1260', 'SiIII_1206', 'CII_1036', 'CIII_977', 'OVI_1032', 'OVI_1032', 'HI_1215', 'HI_1215', 'HI_1215']

binning=1.23

for i in range(len(lines)):

    print(f'{lines[i]} : {cos_siglevel(eq_w(lines[i],line_name[i]),cen_wave[i],b[i],snpix[i],binning=binning)}')


os.remove('temp_param_file.txt')

quit()






W=eq_w('OVI_1032_1','OVI_1032')
binning=1.23

sl=cos_siglevel(W=l[5],wavelength=l[0],b=l[7],snpix=l[4],binning=binning)



line_param=[[1390.60,  "OVI 1032",   0.347567, 21.5, 18.0,  155,   9,  30.2],
            [1390.92,  "OVI 1032",   0.347878, 14.7, 19.0,   82,   7,  19.0],
            [1396.88,  "CII 1036" ,  0.347903,  6.2, 15.0,   31,  7,  6.4],
            [1416.42,  "Lya 1215",   0.165152, 74.5, 20.0,  686,  40,  45.1],
            [1467.32,  "OVI 1032",   0.421923, 24.1, 21.0,  151,   7,  25.8],
            [1534.87,  "SiIV 1402",  0.094169,  3.4, 21.0,   33,  20,  35.5],
            [1626.31,  "SiIII 1206", 0.347982,  8.4, 11.0,   76,  13,   9.3],
            [1690.40,  "NV 1238",    0.364515,  5.3, 12.0,   69,  23,  30.4],
            [1698.95,  "SiII 1260",  0.347921,  6.6, 12.0,   49,  22,   5.7]
            ]

# WAVE    LINE           Z      sig  S/N   ew  ewr   b

# print(cos_siglevel(eq_W1,wavelength1,b1,snpix=snpix1,binning=binning1))
# print(cos_siglevel(eq_W2,wavelength2,b2,snpix=snpix2,binning=binning1))

binning=1.23

for l in line_param:

    sl=cos_siglevel(W=l[5],wavelength=l[0],b=l[7],snpix=l[4],binning=binning)
    print(f'{sl:.2f} : {l[3]} : {sl/l[3]:.2f}')


