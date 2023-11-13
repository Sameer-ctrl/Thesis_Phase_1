import matplotlib.pyplot as plt
from numpy import loadtxt,ceil, unique
from astropy.io import ascii
import os

plt.style.use('my_style.mpl')

qso='sbs1108'
z_abs=0.463207
vlim=300
n_col=4
vsep=0

spec=ascii.read(f'Data/IGM_Danforth_Data/Cont_norm_spectra/{qso}_cont_norm.asc')
wave_spec=spec['WAVE']
flux_spec=spec['FLUX']
err_spec=spec['ERROR']

qso_list=loadtxt('Data/qso_list.txt',dtype=str)
qso_dict=dict(zip(qso_list[:,1],qso_list[:,0]))
qso_label=qso_dict[qso]


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


# lines=[]

# for file in files:

#     if file[:2]!='Ly':

#         a=file.split('_')

#         if len(a)==3:
#             lines.append(f'{a[0]}_{a[1]}')
        
#         else:
#             b=a[1].split('.')[0]
#             lines.append(f'{a[0]}_{b}')
    
#     else:
#         a=file.split('_')

#         if len(a)==2:
#             lines.append(a[0])
        
#         else:
#             lines.append(a[0].split('.')[0])


# lines=unique(lines)

# print(lines)
# quit()

lines=['CIII_977', 'CII_1036', 'HI_1215', 'HI_1025', 'HI_972', 'HI_937', 'HI_949' , 'HI_937', 'HI_930', 'HI_926', 'HI_923', 'HI_920'
 ,'NIII_989', 'OI_988', 'OVI_1032', 'OVI_1038', 'SiIII_1206', 'SiII_1190', 'SiII_1193']

# lines=['HI_1215','OVI_1032','SiIII_1206','HI_1025','OVI_1038','CIII_977','HI_972','HI_949','HI_937','HI_930','HI_926','HI_923']
# lines=['HI_1215','OVI_1032','CIII_977','HI_1025','OVI_1038','SiIII_1206','HI_972','HI_949','HI_937']
# line_label={'HI_1215':['H','I','1215'],'HI_1025':['H','I','1025'],'HI_972':['H','I','972'],'HI_949':['H','I','949'],'HI_937':['H','I','937'],'HI_930':['H','I','930'],'OVI_1032':['O','VI','1032'],'CII_1036':['C','II','1036'],'CII_1334':['C','II','1334'],'HI_1025':['H','I','1025'],'OVI_1038':['O','VI','1038'],'OI_1302':['O','I','1302'],'CIII_977':['C','III','977'],'CIV_1548':['C','IV','1548'],'CIV_1550':['C','IV','1550'],'HI_972':['H','I','972'],'SiII_1260':['Si','II','1260'],'SiII_1190':['Si','II','1190'],'SiII_1193':['Si','II','1193'],'SiII_1304':['Si','II','1304'],'SiII_1526':['Si','II','1526'],'SiIII_1206':['Si','III','1206'],'SiIV_1393':['Si','IV','1393'],'SiIV_1402':['Si','IV','1402'],'NII_1083':['N','II','1083'],'NIII_989':['N','III','989'],'NV_1238':['N','V','1238'],'NV_1242':['N','V','1242'],'PII_1152':['P','II','1152']}
# lines=['HI_1215','HI_1025','HI_972','HI_949','HI_937','HI_930','HI_926','HI_923','HI_920','HI_919','HI_918','HI_917','HI_916']

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


def abs_line_plot(line):

    data=loadtxt(f'{file_path}/{line}.txt',comments='!')

    cen_wave_rest=rest_wave[line]
    cen_wave_obs=cen_wave_rest*(1+z_abs)

    wave=data[:,0]
    cont1=data[:,3]
    v1=3*(10**5)*((wave**2-(cen_wave_obs**2))/(wave**2+(cen_wave_obs**2)))
    v2=3*(10**5)*((wave_spec**2-(cen_wave_obs**2))/(wave_spec**2+(cen_wave_obs**2)))

    plt.step(v2,flux_spec,c='green',lw=2,label=r'$\mathbf{Flux}$')

    vpfit_chunks=file_group(line)

    n=len(vpfit_chunks)-1
    m=0

    for f in vpfit_chunks:
        splitted=f.split('_')
        if len(splitted)>2:
            if splitted[2][:4]=='cont':
                m+=1

    if n-m>0:

        for i in range(n-m):

            comp=i+1
            data=loadtxt(f'{file_path}/{line}_{comp}.txt',comments='!')

            wave=data[:,0]
            cont=data[:,3]

            v=3*(10**5)*((wave**2-(cen_wave_obs**2))/(wave**2+(cen_wave_obs**2)))

            plt.plot(v,cont,ls='--')

    if m>0:

        for i in range(m):

            if m==1:
                data=loadtxt(f'{file_path}/{line}_cont.txt',comments='!')

            else:
                data=loadtxt(f'{file_path}/{line}_cont{i+1}.txt',comments='!')

            wave=data[:,0]
            cont=data[:,3]

            v=3*(10**5)*((wave**2-(cen_wave_obs**2))/(wave**2+(cen_wave_obs**2)))

            plt.plot(v,cont,ls='--',c='black')



    plt.hlines(1,-5000,5000,ls='--',lw=1,color='black')
    plt.hlines(0,-5000,5000,ls='--',lw=1,color='black')
    plt.vlines(0,-1,2,ls='--',lw=1,color='black')

    plt.plot(v1,cont1,c='red',label=r'$\mathbf{Voigt \ profile \ fit}$',lw=3)
    plt.step(v2,err_spec,c='#ffb3ff',label=r'$\mathbf{Error}$',lw=2)
    
    for i in range(len(line)):
        if line[i]=='_':
            break

    plt.xlim(-vlim,vlim)
    plt.ylim(-0.1,1.6)
    line_name=line_label(line)
    plt.text(-260,0.21,f'{{\\fontsize{{25pt}}{{3em}}\selectfont{{}}$\mathbf{{{line_name[0]}}}$}} {{\\fontsize{{17pt}}{{3em}}\selectfont{{}}$\mathbf{{{line_name[1]}}}$}} {{\\fontsize{{25pt}}{{3em}}\selectfont{{}}$\mathbf{{{line_name[2]}}}$}}')
    plt.yticks([0,0.5,1,1.5])
    plt.xticks([-300,-200,-100,0,100,200,300])
    

fig=plt.figure(figsize=(40,20),dpi=300)

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
plt.legend(bbox_to_anchor=(0.51,1.03),bbox_transform=plt.gcf().transFigure, loc='center',ncols=3,fontsize=30)
plt.text(0.38, 1.08, f'$\mathbf{{{qso_label} \ (z_{{abs}}={z_abs})}}$', fontsize=40, transform=plt.gcf().transFigure)
plt.savefig(f'Files_n_figures/{qso_label}_z={z_abs}_sys_plot.png')
plt.savefig(f'../VPfit/{qso}/z={z_abs}/{qso_label}_z={z_abs}_sys_plot.png')

# plt.show()  