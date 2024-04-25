from  numpy import *
from scipy.constants import G, m_e,m_p,speed_of_light,parsec
import os
from astropy.io import ascii,fits
import matplotlib.pyplot as plt
from astropy.table import Table,Column

lyman_alpha={'3c263': 1386.783, 'pks0637': 1723.26, 'pg1424': 1394.5, 'pg0003': 1728.589, 'pg1216': 1558.837, 's135712': 1334.647, '1es1553': 1443.929, 'sbs1108': 1778.777, 'pg1222': 1675.666, 'pg1116': 1384.073, 'h1821': 1489.173, 'pg1121': 1449.557, 'pks0405': 1418.839}

selected_plot = None
selected_points = []

selected_x = []
selected_y = []

x_min, x_max, y_min, y_max = None, None, None, None

i=0

def write_exclude_wave(qso):

    data=ascii.read(f'Data/IGM_Danforth_Data/Cont_norm_spectra/{qso}_cont_norm.asc')
    wave=data['WAVE']
    flux=data['FLUX']

    # selected_points = []
    # selected_plot = None

    # selected_x = []
    # selected_y = []

    # x_min, x_max, y_min, y_max = None, None, None, None

    # i=0

    def onclick(event): 

        global i, selected_points

        if event.key == 'a': 
        
            i+=1
            selected_points.append((event.xdata, event.ydata))
            print(f'{i} - A : {event.xdata}, {event.ydata}')
            ax.scatter(event.xdata, event.ydata, color='red')  # Highlight selected point

            fig.canvas.draw()  # Create a figure and axis

        if event.key=='d':
            i-=1
            print(f'{i} - D : {selected_points[-1][0]}, {selected_points[-1][1]}')
            selected_points.pop()
            
        update_plot(event.key)


    def update_plot(key=0):

        global selected_plot, x_min, x_max, y_min, y_max

        if key=='q':
            plt.show(block=False)

        else:

            if plt.gca().has_data():
                x_min, x_max = plt.xlim()
                y_min, y_max = plt.ylim()

            plt.clf()  # Clear the current figure

            plt.step(wave,flux)
            plt.vlines(lyman_alpha[qso],-1,5,ls='--',lw=3)

            if selected_points:
                
                x_values, y_values = zip(*selected_points)

                selected_plot = plt.scatter(x_values, y_values, color='red')

                if len(selected_points)>1:

                    for k in range(len(selected_points)-1):
                        if k%2==0:
                            plt.hlines(y_values[k],x_values[k],x_values[k+1],color='green')

            if x_min and x_max and y_min and y_max:
                plt.xlim(x_min, x_max)
                plt.ylim(y_min, y_max)

            plt.draw()


    fig,ax=plt.subplots()

    ax.set_xlabel('Wavelength')
    ax.set_ylabel('Flux')
    plt.step(wave,flux)
    plt.vlines(lyman_alpha[qso],-1,5,ls='--',lw=3)
    cid=fig.canvas.mpl_connect('key_press_event', onclick)

    plt.show()


    x=array([p[0] for p in selected_points])

    wave_exc_l=zeros(int(len(x)/2))
    wave_exc_r=zeros(int(len(x)/2))
    l=0
    r=0

    for j in range(len(x)):

        if j%2==0:
            wave_exc_l[l]=round(x[j],3)
            l+=1
        
        else:
            wave_exc_r[r]=round(x[j],3)
            r+=1

    tab=Table()

    wave_exc_l_col=Column(name='WAVE1',data=wave_exc_l)
    wave_exc_r_col=Column(name='WAVE2',data=wave_exc_r)

    tab.add_columns([wave_exc_l_col,wave_exc_r_col])
    tab.write(f'Data/IGM_Danforth_Data/Excluded_wavelengths/{qso}_excluded_wavelength.asc', format='ascii', overwrite=True)


def redshift_path(qso,wave_min=1220,v_lim=5000):

    file_systems=open(f'Data/IGM_Danforth_Data/Systems/{qso}_igm-systems.txt','r')
    z_em=float(file_systems.readlines()[16].split(' ')[1])

    data=ascii.read(f'Data/IGM_Danforth_Data/Cont_norm_spectra/{qso}_cont_norm.asc')
    wave=data['WAVE']

    excluded_wave=ascii.read(f'Data/IGM_Danforth_Data/Excluded_wavelengths/{qso}_excluded_wavelength.asc')

    wave_l=excluded_wave['WAVE1']
    wave_r=excluded_wave['WAVE2']

    v_z=lambda z : 3e5*(((1+z)**2-1)/((1+z)**2+1))  # v at z
    z_v=lambda v : sqrt((1+((v)/3e5))/(1-((v)/3e5)))-1      # z at v

    z_lim=z_v(v_z(z_em)-v_lim)

    rest_wave=1215.6701
    wave_max=min([(1+z_lim)*rest_wave,wave[-1]])

    dzb=round(sum((wave_r-wave_l)/rest_wave),3)

    dz=round((wave_max-wave_min)/rest_wave,3)
    dz_unblocked=round(dz-dzb,3)

    print(dz,dzb,dz_unblocked)


def plot_excluded_region(qso,y=1.25,dy=0.1):

    data=ascii.read(f'Data/IGM_Danforth_Data/Cont_norm_spectra/{qso}_cont_norm.asc')
    wave=data['WAVE']
    flux=data['FLUX']

    excluded_wave=ascii.read(f'Data/IGM_Danforth_Data/Excluded_wavelengths/{qso}_excluded_wavelength.asc')

    wave_l=excluded_wave['WAVE1']
    wave_r=excluded_wave['WAVE2']

    plt.step(wave,flux)

    for i in range(len(wave_l)):
        plt.hlines(y,wave_l[i],wave_r[i],color='red',lw=2)
        plt.vlines(wave_l[i],y-dy,y+dy,color='red',lw=2)
        plt.vlines(wave_r[i],y-dy,y+dy,color='red',lw=2)

    plt.xlabel('Wavelength')
    plt.ylabel('Flux')
    plt.ylim(-0.3,3)
    plt.show()
        

qso='pg1121'


# pks0405 : *
# pg1116  : *
# sbs1108 : *
# pg1121  : *
# s135712 : 5
# pg1424  : 4
# pg1222  : 4
# pg0003  : 4,3
# pks0637 : 3,3
# h1821   : 3,3
# 1es1553 : 3
# 3c263   : 3
# pg1216  : 3


write_exclude_wave(qso)
redshift_path(qso)
# plot_excluded_region(qso)
 