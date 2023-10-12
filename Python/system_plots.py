import matplotlib.pyplot as plt
from numpy import loadtxt,ceil, unique
from astropy.io import ascii
import os

plt.style.use('my_style.mpl')

qso='pks0405'
dec=''
z_abs=0.167125 

vlim=350

file_path=f'../VPfit/{qso}/z={z_abs}/VPfit_chunks'
files=os.listdir(f'{file_path}')
data=loadtxt('Data/rest_wave.txt',dtype=str)

ion=data[:,1]
wave=data[:,0].astype(float)

rest_wave={}

for i in range(len(ion)):
    rest_wave.update({ion[i]:wave[i]})


def file_group(x):

    grouped=[]
    n=len(x)

    for i in files:   
        if i[0:n]==x:
            grouped.append(i)
    
    return grouped

lines=[]

for file in files:
    a=file.split('_')
    
    if len(a)==3:
        lines.append(f'{a[0]}_{a[1]}')
    
    else:
        b=a[1].split('.')[0]
        lines.append(f'{a[0]}_{b}')

lines=unique(lines)


# lines=['HI_1215','OVI_1032','CII_1036','CII_1334','HI_1025','OVI_1038','SiII_1260','SiIII_1206','SiII_1193','SiII_1190','SiIV_1393','NII_1083','NIII_989','NV_1238','NV_1242']
# line_label={'HI_1215':['H','I','1215'],'HI_1025':['H','I','1025'],'HI_972':['H','I','972'],'HI_949':['H','I','949'],'HI_937':['H','I','937'],'HI_930':['H','I','930'],'OVI_1032':['O','VI','1032'],'CII_1036':['C','II','1036'],'CII_1334':['C','II','1334'],'HI_1025':['H','I','1025'],'OVI_1038':['O','VI','1038'],'OI_1302':['O','I','1302'],'CIII_977':['C','III','977'],'CIV_1548':['C','IV','1548'],'CIV_1550':['C','IV','1550'],'HI_972':['H','I','972'],'SiII_1260':['Si','II','1260'],'SiII_1190':['Si','II','1190'],'SiII_1193':['Si','II','1193'],'SiII_1304':['Si','II','1304'],'SiII_1526':['Si','II','1526'],'SiIII_1206':['Si','III','1206'],'SiIV_1393':['Si','IV','1393'],'SiIV_1402':['Si','IV','1402'],'NII_1083':['N','II','1083'],'NIII_989':['N','III','989'],'NV_1238':['N','V','1238'],'NV_1242':['N','V','1242'],'PII_1152':['P','II','1152']}
# lines=['HI_1215','HI_1025','HI_972','HI_949','HI_937','HI_930','HI_926','HI_923','HI_920','HI_919','HI_918']

def line_label(line):

    for i,s in enumerate(line):
        if s.isupper():
            if i==0:
                pass
            else:
                break

    ion=line[:i]
    n,wave=line[i:].split('_')

    return ion,n,wave


n=len(lines)

# f,ax_sub=plt.subplots(n,1,figsize=(15,45))

# quit()


def abs_line_plot(line):

    data=loadtxt(f'{file_path}/{line}.txt',comments='!')
    spec=ascii.read(f'Data/IGM_Danforth_Data/Cont_norm_spectra/{qso}_cont_norm.asc')
    wave_spec=spec['WAVE']
    flux_spec=spec['FLUX']
    err_spec=spec['ERROR']

    cen_wave_rest=rest_wave[line]
    cen_wave_obs=cen_wave_rest*(1+z_abs)

    wave=data[:,0]
    # flux=data[:,1]
    # error=data[:,2]
    cont1=data[:,3]
    v1=3*(10**5)*((wave**2-(cen_wave_obs**2))/(wave**2+(cen_wave_obs**2)))
    v2=3*(10**5)*((wave_spec**2-(cen_wave_obs**2))/(wave_spec**2+(cen_wave_obs**2)))

    # plt.step(v1,flux,c='black')
    plt.step(v2,flux_spec,c='green',lw=2,label=r'$\mathbf{Flux}$')

    n=len(file_group(line))-1
    if n>0:
        for i in range(n):

            comp=i+1
            data=loadtxt(f'{file_path}/{line}_{comp}.txt',comments='!')

            wave=data[:,0]
            cont=data[:,3]
            v=3*(10**5)*((wave**2-(cen_wave_obs**2))/(wave**2+(cen_wave_obs**2)))

            # plt.plot(v,cont,label=f'{comp}',ls='--')
            plt.plot(v,cont,ls='--')

    plt.hlines(1,-5000,5000,ls='--',lw=1,color='black')
    plt.hlines(0,-5000,5000,ls='--',lw=1,color='black')
    plt.vlines(0,-0.1,1.7,ls='--',lw=1,color='black')
    plt.plot(v1,cont1,c='red',label=r'$\mathbf{Voigt \ profile \ fit}$',lw=3)
    plt.step(v2,err_spec,c='#ffb3ff',label=r'$\mathbf{Error}$',lw=2)
    
    for i in range(len(line)):
        if line[i]=='_':
            break

    plt.xlim(-vlim,vlim)
    plt.ylim(-0.1,1.7)
    line_name=line_label(line)
    plt.text(-260,0.21,f'{{\\fontsize{{25pt}}{{3em}}\selectfont{{}}$\mathbf{{{line_name[0]}}}$}} {{\\fontsize{{17pt}}{{3em}}\selectfont{{}}$\mathbf{{{line_name[1]}}}$}} {{\\fontsize{{25pt}}{{3em}}\selectfont{{}}$\mathbf{{{line_name[2]}}}$}}')
    plt.yticks([0,0.5,1,1.5])
    

fig=plt.figure(figsize=(30,18),dpi=300)

n_col=5

ax1=plt.subplot(int(ceil(n/n_col)),n_col,1)

for i,line in enumerate(lines):

    abs_line_plot(line)

    if i%n_col!=0:
        plt.tick_params('y', labelleft=False)

    else:
        plt.tick_params('y', labelsize=25)

    if i < n-n_col:
        plt.tick_params('x', labelbottom=False)
    
    else:
        plt.tick_params('x', labelbottom=True,labelsize=25)

    if i!=n-1:
        ax=plt.subplot(int(ceil(n/n_col)),n_col,i+2,sharex=ax1, sharey=ax1)


fig.supxlabel(r'$\mathbf{V} \ \mathbf{(km \ \ s^{-1})}$',fontsize=30,y=-0.02)  #y=0.18  (y=0 for lyman)
fig.supylabel(r'$\mathbf{Continuum \ Normalized \ Flux} $',fontsize=30,x=0.08, y=0.52) #x=0.05, y=0.62 (x=0.05, y=0.55 for lyman)
plt.subplots_adjust(hspace=0,top=0.99,bottom=0.07,wspace=0)
plt.legend(bbox_to_anchor=(-1.5,4.4), loc='upper center',ncols=3,fontsize=30)
plt.text(0.4, 1.1, f'$\mathbf{{{qso.upper()} \ (z_{{abs}}={z_abs})}}$', fontsize=40, transform=plt.gcf().transFigure)
plt.savefig(f'Files_n_figures/{qso}_sys_plot.png')
plt.savefig(f'../VPfit/{qso}/{qso}_sys_plot.png')

# plt.show()  



# quit()


# t = np.arange(0.01, 5.0, 0.01)
# s1 = np.sin(2 * np.pi * t)
# s2 = np.exp(-t)
# s3 = np.sin(4 * np.pi * t)

# ax1 = plt.subplot(311)
# plt.plot(t, s1)
# plt.tick_params('x', labelsize=6)

# # share x only
# ax2 = plt.subplot(312, sharex=ax1, sharey=ax1)
# plt.plot(t, s2)
# # make these tick labels invisible
# plt.tick_params('x', labelbottom=False)

# # share x and y
# ax3 = plt.subplot(313, sharex=ax1, sharey=ax1)
# plt.plot(t, s3)
# plt.xlim(0.01, 5.0)
# plt.show()