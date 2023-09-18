import matplotlib.pyplot as plt
from numpy import loadtxt,ceil
from astropy.io import ascii
import os

plt.style.use('my_style.mpl')

file_path='Data/VPfit_fits_rebinned/Metals_HI'

files=os.listdir(f'{file_path}')
z_abs=0.347579

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

wave_obs={rest_wave['HI_1215']:[1638.23436,1637.23542,1638.61559],rest_wave['HI_1025']:[1382.36100,1382.55866,1382.79019],rest_wave['HI_972']:[1310.56725],rest_wave['OVI_1032']:[1390.60515,1390.91597],rest_wave['OVI_1038']:[1398.27370,1398.58624],rest_wave['CII_1036']:[1396.94145],rest_wave['CIII_977']:[1316.52808,1316.94620],rest_wave['SiII_1260']:[1698.99935],rest_wave['SiIII_1206']:[1626.31495]}
# wave_rest=[1215,1025,972,1032,1038,1036,977,1260,1206]

# z=[]

# for i in wave_obs:
#     for j in wave_obs[i]:
#         z.append((j-i)/i)

# print(z)

# plt.hist(z,bins='auto')
# plt.show()
lines=['HI_1215','OVI_1032','CII_1036','HI_1025','OVI_1038','CIII_977','HI_972','SiII_1260','SiIII_1206']
line_label={'HI_1215':['H','I','1215'],'OVI_1032':['O','VI','1032'],'CII_1036':['C','II','1036'],'HI_1025':['H','I','1025'],'OVI_1038':['O','VI','1038'],'CIII_977':['C','III','977'],'HI_972':['H','I','972'],'SiII_1260':['Si','II','1260'],'SiIII_1206':['Si','III','1206']}
# lines=['HI_1215','HI_1025','HI_972','HI_949','HI_937','HI_930','HI_926','HI_923','HI_920','HI_919','HI_918']

n=len(lines)

# f,ax_sub=plt.subplots(n,1,figsize=(15,45))

# quit()


class abs_line():

    def __init__(self,line):

        data=loadtxt(f'{file_path}/{line}.txt',comments='!')
        spec=ascii.read('Data/PG0003+158_rebinned_cont_norm.asc')
        wave_spec=spec['WAVE']
        flux_spec=spec['FLUX']
        err_spec=spec['ERROR']

        cen_wave_rest=rest_wave[line]
        cen_wave_obs=cen_wave_rest*(1+z_abs)

        wave=data[:,0]
        # flux=data[:,1]
        error=data[:,2]
        cont1=data[:,3]
        v1=3*(10**5)*((wave**2-(cen_wave_obs**2))/(wave**2+(cen_wave_obs**2)))
        v2=3*(10**5)*((wave_spec**2-(cen_wave_obs**2))/(wave_spec**2+(cen_wave_obs**2)))

        # plt.step(v1,flux,c='black')
        plt.step(v2,flux_spec,c='green',lw=2,label=r'$\mathbf{Flux}$')

        n=len(file_group(line))-1
        # component=['I','II','III','IV']
        if n>0:
            for i in range(n):

                comp=i+1
                data=loadtxt(f'{file_path}/{line}_{comp}.txt',comments='!')

                wave=data[:,0]
                cont=data[:,3]
                v=3*(10**5)*((wave**2-(cen_wave_obs**2))/(wave**2+(cen_wave_obs**2)))

                plt.plot(v,cont,label=f'{comp}',ls='--')

        plt.hlines(1,-350,350,ls='--',lw=1,color='black')
        plt.hlines(0,-350,350,ls='--',lw=1,color='black')
        plt.vlines(0,-0.1,1.7,ls='--',lw=1,color='black')
        plt.plot(v1,cont1,c='red',label=r'$\mathbf{Voigt \ profile \ fit}$',lw=3)
        plt.step(v2,err_spec,c='#ffb3ff',label=r'$\mathbf{Error}$',lw=2)
        
        for i in range(len(line)):
            if line[i]=='_':
                break

        # plt.ylabel(f'{line[:i]} {line[i+1:]}',fontsize=20)
        plt.xlim(-350,350)
        plt.ylim(-0.1,1.7)
        line_name=line_label[line]
        # plt.text(-260,0.21,r'{\fontsize{25pt}{3em}\selectfont{}$\mathbf{1206}$} {\fontsize{17pt}{3em}\selectfont{}$\mathbf{VI}$')
        plt.text(-260,0.21,f'{{\\fontsize{{25pt}}{{3em}}\selectfont{{}}$\mathbf{{{line_name[0]}}}$}} {{\\fontsize{{17pt}}{{3em}}\selectfont{{}}$\mathbf{{{line_name[1]}}}$}} {{\\fontsize{{25pt}}{{3em}}\selectfont{{}}$\mathbf{{{line_name[2]}}}$}}')
        # plt.text(-260,0.11,f'$\mathbf{{{line[:i]} \ \ {line[i+1:]}}}$',fontsize=20)
        # plt.text(-250,0.3,r'$\textbf{}$'
        # plt.yticks(fontsize=15)
        # plt.tick_params(axis="x", labelsize=15)
        # plt.tick_params(axis="y", labelsize=15)
    

fig=plt.figure(figsize=(24,12),dpi=300)
# ax = fig.add_subplot(111,frameon=False)
# ax.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)



ax1=plt.subplot(int(ceil(n/3)),3,1)
# k=0

for i in range(n):   

    a=abs_line(lines[i])
    plt.tick_params('x', labelbottom=False)
    plt.yticks([0,0.5,1,1.5])
    # plt.xticks([-300,-200,-100,0,100,200,300])

    if i!=n-1:

        if i%3!=0:
            plt.tick_params('y', labelleft=False)
        
        else:
            plt.tick_params('y', labelsize=25)
        
        if i==n-2 or i==n-3:
           plt.tick_params('x', labelbottom=True,labelsize=25) 
        
        ax=plt.subplot(int(ceil(n/3)),3,i+2,sharex=ax1, sharey=ax1)
        
    else:
        plt.tick_params('x', labelbottom=True,labelsize=25)
        if n%3==0:
            plt.tick_params('y', labelleft=False)
        

fig.supxlabel(r'$\mathbf{V} \ \mathbf{(km \ \ s^{-1})}$',fontsize=30,y=-0.02)  #y=0.18  (y=0 for lyman)
fig.supylabel(r'$\mathbf{Continuum \ Normalized \ Flux} $',fontsize=30,x=0.08, y=0.52) #x=0.05, y=0.62 (x=0.05, y=0.55 for lyman)
plt.subplots_adjust(hspace=0,top=0.99,bottom=0.07,wspace=0)
plt.legend(bbox_to_anchor=(-0.5,3.38), loc='upper center',ncols=3,fontsize=30)
# plt.legend(bbox_to_anchor=(1,6.5), loc='upper center',ncols=3)  #(1,6.5) for lyman 
plt.savefig('Files_n_figures/sys_plot_rebinned_3.png')
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