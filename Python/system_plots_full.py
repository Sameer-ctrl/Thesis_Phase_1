import matplotlib.pyplot as plt
from numpy import *
from scipy.integrate import simpson
from astropy.io import ascii
import os

plt.style.use('my_style.mpl')


qso='3c57'
z_abs=0.077430
vlim=350
n_col=6
lw1=1.5
vlim_aod=[-100,100]

spec=ascii.read(f'Data/IGM_Danforth_Data/Cont_norm_spectra/{qso}_cont_norm.asc')
wave_spec=spec['WAVE']
flux_spec=spec['FLUX']
err_spec=spec['ERROR']

qso_list=loadtxt('Data/qso_list.txt',dtype=str)
qso_dict=dict(zip(qso_list[:,1],qso_list[:,0]))
qso_label=qso_dict[qso]


file_path=f'../VPfit/{qso}/z={z_abs:.6f}/VPfit_chunks'
files=os.listdir(f'{file_path}')
data=loadtxt('Data/lines_system_plots.txt',dtype=str)

line_atom=data[:,0]
wave_atom=data[:,1].astype(float)
f_atom=data[:,2].astype(float)

rest_wave={}
f={}

for i in range(len(line_atom)):
    rest_wave.update({line_atom[i]:wave_atom[i]})
    f.update({line_atom[i]:f_atom[i]})


def file_group(x):

    grouped=[]
    n=len(x)

    for i in files:   
        if i[0:n]==x:
            grouped.append(i)
    
    return grouped


def lines_all():

    lines=[]

    for file in files:

        if file[:2]!='Ly':

            a=file.split('_')

            if len(a)==3:
                lines.append(f'{a[0]}_{a[1]}')
            
            else:
                b=a[1].split('.')[0]
                lines.append(f'{a[0]}_{b}')
        
        else:
            a=file.split('_')

            if len(a)==2:
                lines.append(a[0])
            
            else:
                lines.append(a[0].split('.')[0])


    lines=unique(lines)

    return lines


def tick_pos(file):

    with open(file) as f:
        data=f.readlines()

    i=0 
    tick_wave=[]

    for line in data:
        
        if line[0]=='!':
            i+=1

            if i>2:
                tick_wave.append(float(line[:-1].split()[1]))

    return array(tick_wave)


# lines=lines_all()
lines=line_atom
# print(lines)
# quit()
# lines=['HI_1215', 'OVI_1032','CIII_977', 'HI_1025' , 'OVI_1038', 'CII_1036', 'HI_972' ,'SiIII_1206', 'SiII_1260']
# lines=['HI_1215', 'HI_1025', 'HI_972', 'OVI_1032', 'OVI_1038','CIII_977' , 'CII_1036' ,'SiIII_1206', 'SiII_1260']
n=len(lines)
line_vshift=dict(zip(lines,zeros(n)))
# line_vshift['HI_1215']=-10
# line_vshift['HI_972']=4
# line_vshift['OVI_1038']=3
# line_vshift['SiIII_1206']=2 
# line_vshift['OVI_1038']=-5
# line_vshift['HI_949']=-3
# line_vshift['HI_937']=-12
# line_vshift['HI_930']=-4
# line_vshift['SiII_1193']=-4
# line_vshift['CIV_1548']=-5
# line_vshift['CIV_1550']=9
# line_vshift['SiIV_1402']=-5
# line_vshift['SiIV_1393']=-5
# line_vshift['SiII_1193']=-6
# line_vshift['SiII_1190']=-5
# line_vshift['HI_972']=-10
# line_vshift['HI_1215']=-3
# line_vshift['HI_916']=-4
# line_vshift['HI_917']=-3
# line_vshift['HI_918']=-2
# line_vshift['HI_972']=-1
# line_vshift['HI_949']=-1
# line_vshift['HI_937']=-3
# line_vshift['HI_930']=-4
# line_vshift['HI_926']=-4
# line_vshift['HI_923']=-2
# line_vshift['HI_920']=-2
# line_vshift['Ly14']=-1
# line_vshift['Ly15']=-3
# line_vshift['Ly16']=-3
# line_vshift['SiII_1193']=-3
# line_vshift['SiIII_1206']=-5


# lines=['HI_1215','OVI_1032', 'CIII_977','HI_1025','OVI_1038', 'CII_1036',  'HI_972'  
#  ,'SiIII_1206', 'SiII_1260']

# quit()

def line_label(line):

    if line[:2]=='Ly':
        ion=line
        n=''
        wave=''
    
    else:

        for i,s in enumerate(line):
            if s.isupper():
                if i==0:
                    pass
                else:
                    break

        ion=line[:i]
        n,wave=line[i:].split('_')

    return ion,n,wave


def abs_line_plot(line):

    cen_wave_rest=rest_wave[line]
    cen_wave_obs=cen_wave_rest*(1+z_abs)
    vshift=line_vshift[line]

    print(f'{line}  :  {cen_wave_obs}')

    if os.path.exists(f'{file_path}/{line}.txt'):

        data=loadtxt(f'{file_path}/{line}.txt',comments='!')


        wave=data[:,0]
        cont1=data[:,3]
    
        v1=3*(10**5)*((wave**2-(cen_wave_obs**2))/(wave**2+(cen_wave_obs**2)))
        v2=3*(10**5)*((wave_spec**2-(cen_wave_obs**2))/(wave_spec**2+(cen_wave_obs**2)))

        plt.step(v2-vshift,flux_spec,c='green',lw=lw1+2,label=r'$\mathbf{Flux}$')

        vpfit_chunks=file_group(line)

        n=len(vpfit_chunks)-1
        m=0

        for file in vpfit_chunks:
            splitted=file.split('_')

            if file[:2]=='Ly':
                if len(splitted)>1:
                    if splitted[1][:4]=='cont':
                        m+=1
                
            else:

                if len(splitted)>2:
                    if splitted[2][:4]=='cont':
                        m+=1
            
        
        if n-m>0:

            for i in range(n-m):

                comp=i+1
                data=loadtxt(f'{file_path}/{line}_{comp}.txt',comments='!')
                tick_wave=tick_pos(f'{file_path}/{line}_{comp}.txt')

                wave=data[:,0]
                cont=data[:,3]

                v=3*(10**5)*((wave**2-(cen_wave_obs**2))/(wave**2+(cen_wave_obs**2)))
                v_tick=3*(10**5)*((tick_wave**2-(cen_wave_obs**2))/(tick_wave**2+(cen_wave_obs**2)))

                plt.plot(v,cont,ls='--',lw=lw1+2,c='blue')
                plt.vlines(v_tick,1.1,1.2,color='blue',lw=2.5)

        if m>0:

            for i in range(m):

                if m==1:
                    data=loadtxt(f'{file_path}/{line}_cont.txt',comments='!')
                    tick_wave=tick_pos(f'{file_path}/{line}_cont.txt')


                else:
                    data=loadtxt(f'{file_path}/{line}_cont{i+1}.txt',comments='!')
                    tick_wave=tick_pos(f'{file_path}/{line}_cont{i+1}.txt')

                wave=data[:,0]
                cont=data[:,3]

                v=3*(10**5)*((wave**2-(cen_wave_obs**2))/(wave**2+(cen_wave_obs**2)))
                v_tick=3*(10**5)*((tick_wave**2-(cen_wave_obs**2))/(tick_wave**2+(cen_wave_obs**2)))


                plt.plot(v,cont,ls='--',c='orange',lw=lw1+1.5)
                plt.vlines(v_tick,1.1,1.2,color='orange',lw=2.5)


        plt.hlines(1,-5000,5000,ls='--',lw=1,color='black')
        plt.hlines(0,-5000,5000,ls='--',lw=1,color='black')
        plt.vlines(0,-1,2,ls='--',lw=1,color='black')


        tick_wave=tick_pos(f'{file_path}/{line}.txt')
        v_tick=3*(10**5)*((tick_wave**2-(cen_wave_obs**2))/(tick_wave**2+(cen_wave_obs**2)))

        plt.plot(v1,cont1,c='red',label=r'$\mathbf{Voigt \ profile \ fit}$',lw=lw1+3)
        plt.step(v2-vshift,err_spec,c='#ffb3ff',label=r'$\mathbf{Error}$',lw=lw1+2)

        if n-m==0:
            plt.vlines(v_tick,1.1,1.2,color='blue',lw=2.5)


        
        for i in range(len(line)):
            if line[i]=='_':
                break

        

        plt.xlim(-vlim,vlim)
        plt.ylim(-0.1,1.6)
        line_name=line_label(line)
        plt.text(-260,0.21,f'{{\\fontsize{{35pt}}{{3em}}\selectfont{{}}$\mathbf{{{line_name[0]}}}$}} {{\\fontsize{{25pt}}{{3em}}\selectfont{{}}$\mathbf{{{line_name[1]}}}$}} {{\\fontsize{{35pt}}{{3em}}\selectfont{{}}$\mathbf{{{line_name[2]}}}$}}')
        # plt.text(-260,0.21,f'{{\\fontsize{{50pt}}{{3em}}\selectfont{{}}$\mathbf{{{line_name[0]}}}$}} {{\\fontsize{{40pt}}{{3em}}\selectfont{{}}$\mathbf{{{line_name[1]}}}$}} {{\\fontsize{{50pt}}{{3em}}\selectfont{{}}$\mathbf{{{line_name[2]}}}$}}')
        plt.yticks([0,0.5,1,1.5])
        plt.xticks([-300,-200,-100,0,100,200,300])

    else:

        v=3*(10**5)*((wave_spec**2-(cen_wave_obs**2))/(wave_spec**2+(cen_wave_obs**2)))

        plt.step(v-vshift,flux_spec,c='green',lw=lw1+2,label=r'$\mathbf{Flux}$')
        plt.step(v-vshift,err_spec,c='#ffb3ff',label=r'$\mathbf{Error}$',lw=lw1+2)

        plt.hlines(1,-5000,5000,ls='--',lw=1,color='black')
        plt.hlines(0,-5000,5000,ls='--',lw=1,color='black')
        plt.vlines(0,-1,2,ls='--',lw=1,color='black')


        plt.xlim(-vlim,vlim)
        plt.ylim(-0.1,1.6)
        line_name=line_label(line)
        plt.text(-260,0.21,f'{{\\fontsize{{35pt}}{{3em}}\selectfont{{}}$\mathbf{{{line_name[0]}}}$}} {{\\fontsize{{25pt}}{{3em}}\selectfont{{}}$\mathbf{{{line_name[1]}}}$}} {{\\fontsize{{35pt}}{{3em}}\selectfont{{}}$\mathbf{{{line_name[2]}}}$}}')
        plt.yticks([0,0.5,1,1.5])
        plt.xticks([-300,-200,-100,0,100,200,300])


        wave_range_aod=cen_wave_obs*sqrt(((3*(10**5))+array(vlim_aod))/((3*(10**5))-array(vlim_aod)))

        if wave_range_aod[0]>min(wave_spec) and wave_range_aod[1]<max(wave_spec):

            f_line=f[line]
            mask=logical_and(v>vlim_aod[0],v<vlim_aod[1])

            v_aod=v[mask]
            tau_v=-log(flux_spec[mask])

            N=simpson(tau_v,v_aod)/((f_line*cen_wave_rest*2.654E-15))

            col_den_aod=round(log10(abs(N)),1)
            plt.text(180,0.21,f'{{\\fontsize{{20pt}}{{3em}}\selectfont{{}}$\mathbf{{<}}$}} {{\\fontsize{{25pt}}{{3em}}\selectfont{{}}$\mathbf{{{col_den_aod}}}$}}')
            plt.vlines(vlim_aod[0],0.9,1.1,lw=3,color='red')
            plt.vlines(vlim_aod[1],0.9,1.1,lw=3,color='red')
            plt.hlines(1,vlim_aod[0],vlim_aod[1],lw=3,color='red')

        
    

fig=plt.figure(figsize=(40,20),dpi=300)
# fig=plt.figure(figsize=(20,20),dpi=300)

ax1=plt.subplot(int(ceil(n/n_col)),n_col,1)

for i,line in enumerate(lines):

    abs_line_plot(line)

    if i%n_col!=0:
        plt.tick_params('y', labelleft=False)

    else:
        plt.tick_params('y', labelsize=35)

    if i < n-n_col:
        plt.tick_params('x', labelbottom=False)
    
    else:
        plt.tick_params('x', labelbottom=True,labelsize=35)

    if i!=n-1:
        ax=plt.subplot(int(ceil(n/n_col)),n_col,i+2,sharex=ax1, sharey=ax1)


fig.supxlabel(r'$\mathbf{V} \ \mathbf{(km \ \ s^{-1})}$',fontsize=50,y=-0.02)  #y=0.18  (y=0 for lyman)
fig.supylabel(r'$\mathbf{Continuum \ Normalized \ Flux} $',fontsize=50,x=0.08, y=0.52) #x=0.05, y=0.62 (x=0.05, y=0.55 for lyman)
plt.subplots_adjust(hspace=0,top=0.99,bottom=0.07,wspace=0)
plt.legend(bbox_to_anchor=(0.51,1.03),bbox_transform=plt.gcf().transFigure, loc='center',ncols=3,fontsize=30)
plt.text(0.38, 1.08, f'$\mathbf{{{qso_label} \ (z_{{abs}}={z_abs:.6f})}}$', fontsize=40, transform=plt.gcf().transFigure)
plt.savefig(f'Files_n_figures/sys_plots_full/{qso_label}_z={z_abs:.6f}_sys_plot_full.png')
# plt.savefig(f'../VPfit/{qso}/z={z_abs:.6f}/{qso_label}_z={z_abs:.6f}_sys_plot_full.png')
# plt.savefig(f'{qso_label}_z={z_abs:.6f}_sys_plot_full.png')




# plt.show()  