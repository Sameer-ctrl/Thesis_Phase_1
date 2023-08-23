import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.axes import Axes as ax
fig, axs = plt.subplots(nrows = 8, ncols = 2, sharex=True, sharey=True, figsize=(5.7,18))

lamda = np.loadtxt('Lya', usecols = 0, dtype = 'float')
flux  = np.loadtxt('Lya', usecols = 1, dtype = 'float')
error = np.loadtxt('Lya', usecols = 2, dtype = 'float')
cont = np.loadtxt('Lya', usecols = 3, dtype = 'float')
cw= 1215.6701
z_comp1 = 0.1385177451
cwz=(0.1384926549+1)*cw
cws=math.pow(cwz,2)
vel=[]
for i in range(len(lamda)):
    velocity = 300000*(((np.power(lamda[i],2))-cws)/((np.power(lamda[i],2))+cws))
    vel.append(velocity)
f=axs[0,0].plot(vel,flux,color='g',drawstyle='steps-mid', lw=0.7, alpha=0.8)
p=axs[0,0].plot(vel[50:183],cont[50:183],color='red',lw=1.2)
e=axs[0,0].plot(vel,error,color='y', drawstyle='steps-mid',lw=0.8)
axs[0,0].set_xlim(-1.90e+02,1.90e+02)
axs[0,0].set_ylim(-0.1e+00,1.3e+00)
axs[0,0].set_title('HI 1215',fontsize=5, x=0.2, y=0.2 )
axs[0,0].tick_params(axis='both', which='major', labelsize=4.5)
axs[0,0].tick_params(axis='both', which='minor', labelsize=4.5)
axs[0,0].axhline(y=1.00,c="black",linewidth=0.6,ls='dotted')
axs[0,0].axhline(y=-0.1,c="black",linewidth=0.6,ls='dotted')
axs[0,0].axvline(0.5,c="black",linewidth=0.5,ls='dashed')
axs[0,0].spines['bottom'].set_visible(True)
axs[0,0].spines['top'].set_visible(True)
axs[0,0].text(65,0.4,'$λ_{obs}$ = 1384.062', fontsize=5)
   
lamda = np.loadtxt('Lyb', usecols = 0, dtype = 'float')
flux  = np.loadtxt('Lyb', usecols = 1, dtype = 'float')
error = np.loadtxt('Lyb', usecols = 2, dtype = 'float')
cont = np.loadtxt('Lyb', usecols = 3, dtype = 'float')
cw= 1025.7223
z_comp1 = 0.1385177451
cwz=(0.1384926549+1)*cw
cws=math.pow(cwz,2)
vel=[]
for i in range(len(lamda)):
    velocity = 300000*(((np.power(lamda[i],2))-cws)/((np.power(lamda[i],2))+cws))
    vel.append(velocity)
axs[1,0].plot(vel,flux,color='g',drawstyle='steps-mid', lw=0.7, alpha=0.8)
axs[1,0].plot(vel,error,color='y', drawstyle='steps-mid',lw=0.8)
axs[1,0].set_xlim(-1.90e+02,1.90e+02)
axs[1,0].set_ylim(-0.1e+00,1.3e+00)
axs[1,0].set_title('HI 1025',fontsize=5, x=0.2, y=0.2)
axs[1,0].tick_params(axis='both', which='major', labelsize=4.5)
axs[1,0].tick_params(axis='both', which='minor', labelsize=4.5)
axs[1,0].plot(vel[51:132],cont[51:132],color='red',lw=1.2)
axs[1,0].axhline(y=1.00,c="black",linewidth=0.6,ls='dotted')
axs[1,0].axhline(y=-0.1,c="black",linewidth=0.6,ls='dotted')
axs[1,0].axvline(0.5,c="black",linewidth=0.5,ls='dashed')
axs[1,0].spines['bottom'].set_visible(True)
axs[1,0].spines['top'].set_visible(True)
axs[1,0].text(65,0.4,'$λ_{obs}$ = 1167.803', fontsize=5)


lamda = np.loadtxt('CIIa', usecols = 0, dtype = 'float')
flux  = np.loadtxt('CIIa', usecols = 1, dtype = 'float')
error = np.loadtxt('CIIa', usecols = 2, dtype = 'float')
cont = np.loadtxt('CIIa', usecols = 3, dtype = 'float')
cw= 1036.3367
z_comp1 = 0.1384931731
cwz=(0.1384926549+1)*cw
cws=math.pow(cwz,2)
vel=[]
for i in range(len(lamda)):
    velocity = 300000*(((np.power(lamda[i],2))-cws)/((np.power(lamda[i],2))+cws))
    vel.append(velocity)
axs[2,0].plot(vel,flux,color='g',drawstyle='steps-mid', lw=0.7, alpha=0.8)
axs[2,0].plot(vel[51:101],cont[51:101],color='red',lw=1.2)
axs[2,0].plot(vel,error,color='y', drawstyle='steps-mid',lw=0.8)
axs[2,0].set_xlim(-1.90e+02,1.90e+02)
axs[2,0].set_ylim(-0.1e+00,1.3e+00)
axs[2,0].set_title('C II 1036',fontsize=5, x=0.2, y=0.2)
axs[2,0].tick_params(axis='both', which='major', labelsize=4.5)
axs[2,0].tick_params(axis='both', which='minor', labelsize=4.5)
axs[2,0].axhline(y=1.00,c="black",linewidth=0.6,ls='dotted')
axs[2,0].axhline(y=-0.1,c="black",linewidth=0.6,ls='dotted')
axs[2,0].axvline(0.5,c="black",linewidth=0.5,ls='dashed')
axs[2,0].spines['bottom'].set_visible(True)
axs[2,0].spines['top'].set_visible(True)
axs[2,0].text(65,0.4,'$λ_{obs}$ = 1179.862', fontsize=5)

 
lamda = np.loadtxt('CIIb', usecols = 0, dtype = 'float')
flux  = np.loadtxt('CIIb', usecols = 1, dtype = 'float')
error = np.loadtxt('CIIb', usecols = 2, dtype = 'float')
cont = np.loadtxt('CIIb', usecols = 3, dtype = 'float')
cw= 1334.5323
z_comp1 = 0.1384931731
cwz=(0.1384926549+1)*cw
cws=math.pow(cwz,2)
vel=[]
for i in range(len(lamda)):
    velocity = 300000*(((np.power(lamda[i],2))-cws)/((np.power(lamda[i],2))+cws))
    vel.append(velocity)
axs[3,0].plot(vel,flux,color='g',drawstyle='steps-mid', lw=0.7, alpha=0.8)
axs[3,0].plot(vel[50:100],cont[50:100],color='red',lw=1.2)
axs[3,0].plot(vel,error,color='y', drawstyle='steps-mid',lw=0.8)
axs[3,0].set_xlim(-1.90e+02,1.90e+02)
axs[3,0].set_ylim(-0.1e+00,1.3e+00)
axs[3,0].set_title('C II 1334',fontsize=5, x=0.2, y=0.2)
axs[3,0].tick_params(axis='both', which='major', labelsize=4.5)
axs[3,0].tick_params(axis='both', which='minor', labelsize=4.5)
axs[3,0].axhline(y=1.00,c="black",linewidth=0.6,ls='dotted')
axs[3,0].axhline(y=-0.1,c="black",linewidth=0.6,ls='dotted')
axs[3,0].axvline(0.5,c="black",linewidth=0.5,ls='dashed')
axs[3,0].spines['bottom'].set_visible(True)
axs[3,0].spines['top'].set_visible(True)
axs[3,0].text(65,0.4,'$λ_{obs}$ = 1519.356', fontsize=5)

   
lamda = np.loadtxt('CIVa', usecols = 0, dtype = 'float')
flux  = np.loadtxt('CIVa', usecols = 1, dtype = 'float')
error = np.loadtxt('CIVa', usecols = 2, dtype = 'float')
cont = np.loadtxt('CIVa', usecols = 3, dtype = 'float')
cw= 1548.2041
z_comp1 = 0.1384930138
cwz=(0.1384926549+1)*cw
cws=math.pow(cwz,2)
vel=[]
for i in range(len(lamda)):
    velocity = 300000*(((np.power(lamda[i],2))-cws)/((np.power(lamda[i],2))+cws))
    vel.append(velocity)
axs[4,0].plot(vel,flux,color='g',drawstyle='steps-mid', lw=0.7, alpha=0.8)
axs[4,0].plot(vel[51:84],cont[51:84],color='red',lw=1.2)
axs[4,0].plot(vel,error,color='y', drawstyle='steps-mid',lw=0.8)
axs[4,0].set_xlim(-1.90e+02,1.90e+02)
axs[4,0].set_ylim(-0.1e+00,1.3e+00)
axs[4,0].set_title('C IV 1548',fontsize=5, x=0.2, y=0.2)
axs[4,0].tick_params(axis='both', which='major', labelsize=4.5)
axs[4,0].tick_params(axis='both', which='minor', labelsize=4.5)
axs[4,0].axhline(y=1.00,c="black",linewidth=0.6,ls='dotted')
axs[4,0].axhline(y=-0.1,c="black",linewidth=0.6,ls='dotted')
axs[4,0].axvline(0.5,c="black",linewidth=0.5,ls='dashed')
axs[4,0].spines['bottom'].set_visible(True)
axs[4,0].spines['top'].set_visible(True)
axs[4,0].text(65,0.4,'$λ_{obs}$ = 1762.619', fontsize=5)

   
lamda = np.loadtxt('CIVb', usecols = 0, dtype = 'float')
flux  = np.loadtxt('CIVb', usecols = 1, dtype = 'float')
error = np.loadtxt('CIVb', usecols = 2, dtype = 'float')
cont = np.loadtxt('CIVb', usecols = 3, dtype = 'float')
cw= 1550.7812
z_comp1 = 0.1384930138
cwz=(0.1384926549+1)*cw
cws=math.pow(cwz,2)
vel=[]
for i in range(len(lamda)):
    velocity = 300000*(((np.power(lamda[i],2))-cws)/((np.power(lamda[i],2))+cws))
    vel.append(velocity)
axs[5,0].plot(vel,flux,color='g',drawstyle='steps-mid', lw=0.7, alpha=0.8)
axs[5,0].plot(vel[50:94],cont[50:94],color='red',lw=1.2)
axs[5,0].plot(vel,error,color='y', drawstyle='steps-mid',lw=0.8)
axs[5,0].set_xlim(-1.90e+02,1.90e+02)
axs[5,0].set_ylim(-0.1e+00,1.3e+00)
axs[5,0].set_title('C IV 1550',fontsize=5, x=0.2, y=0.2)
axs[5,0].tick_params(axis='both', which='major', labelsize=4.5)
axs[5,0].tick_params(axis='both', which='minor', labelsize=4.5)
axs[5,0].axhline(y=1.00,c="black",linewidth=0.6,ls='dotted')
axs[5,0].axhline(y=-0.1,c="black",linewidth=0.6,ls='dotted')
axs[5,0].axvline(0.5,c="black",linewidth=0.5,ls='dashed')
axs[5,0].spines['bottom'].set_visible(True)
axs[5,0].spines['top'].set_visible(True)
axs[5,0].text(65,0.4,'$λ_{obs}$ = 1765.549', fontsize=5) 
 
   
lamda = np.loadtxt('NII', usecols = 0, dtype = 'float')
flux  = np.loadtxt('NII', usecols = 1, dtype = 'float')
error = np.loadtxt('NII', usecols = 2, dtype = 'float')
cont = np.loadtxt('NII', usecols = 3, dtype = 'float')
cw= 1083.9937
z_comp1 = 0.1385089402
cwz=(0.1384926549+1)*cw
cws=math.pow(cwz,2)
vel=[]
for i in range(len(lamda)):
    velocity = 300000*(((np.power(lamda[i],2))-cws)/((np.power(lamda[i],2))+cws))
    vel.append(velocity)
axs[6,0].plot(vel,flux,color='g',drawstyle='steps-mid', lw=0.7, alpha=0.8)
axs[6,0].plot(vel[51:80],cont[51:80],color='red',lw=1.2)
axs[6,0].plot(vel,error,color='y', drawstyle='steps-mid',lw=0.8)
axs[6,0].set_xlim(-1.90e+02,1.90e+02)
axs[6,0].set_ylim(-0.1e+00,1.3e+00)
axs[6,0].set_title('N II 1083',fontsize=5, x=0.2, y=0.2)
axs[6,0].tick_params(axis='both', which='major', labelsize=4.5)
axs[6,0].tick_params(axis='both', which='minor', labelsize=4.5)
axs[6,0].axhline(y=1.00,c="black",linewidth=0.6,ls='dotted')
axs[6,0].axhline(y=-0.1,c="black",linewidth=0.6,ls='dotted')
axs[6,0].axvline(0.5,c="black",linewidth=0.5,ls='dashed')
axs[6,0].spines['bottom'].set_visible(True)
axs[6,0].spines['top'].set_visible(True)
axs[6,0].text(65,0.4,'$λ_{obs}$ = 1234.136', fontsize=5)

   
   
lamda = np.loadtxt('NVa', usecols = 0, dtype = 'float')
flux  = np.loadtxt('NVa', usecols = 1, dtype = 'float')
error = np.loadtxt('NVa', usecols = 2, dtype = 'float')
cont = np.loadtxt('NVa', usecols = 3, dtype = 'float')
cw= 1238.821
z_comp1 = 0.1385138493
cwz=(0.1384926549+1)*cw
cws=math.pow(cwz,2)
vel=[]
for i in range(len(lamda)):
    velocity = 300000*(((np.power(lamda[i],2))-cws)/((np.power(lamda[i],2))+cws))
    vel.append(velocity)
axs[7,0].plot(vel,flux,color='g',drawstyle='steps-mid', lw=0.7, alpha=0.8)
axs[7,0].plot(vel[51:100],cont[51:100],color='red',lw=1.2)
axs[7,0].plot(vel,error,color='y', drawstyle='steps-mid',lw=0.8)
axs[7,0].set_xlim(-1.900e+02,1.9e+02)
axs[7,0].set_ylim(-0.1e+00,1.3e+00)
axs[7,0].set_title('N V 1238',fontsize=5, x=0.2, y=0.2)
axs[7,0].tick_params(axis='both', which='major', labelsize=4.5)
axs[7,0].tick_params(axis='both', which='minor', labelsize=4.5)
axs[7,0].axhline(y=1.00,c="black",linewidth=0.6,ls='dotted')
axs[7,0].axhline(y=-0.1,c="black",linewidth=0.6,ls='dotted')
axs[7,0].axvline(0.5,c="black",linewidth=0.5,ls='dashed')
axs[7,0].spines['bottom'].set_visible(True)
axs[7,0].spines['top'].set_visible(True)
axs[7,0].text(65,0.4,'$λ_{obs}$ = 1410.414', fontsize=5)

   
lamda = np.loadtxt('NVb', usecols = 0, dtype = 'float')
flux  = np.loadtxt('NVb', usecols = 1, dtype = 'float')
error = np.loadtxt('NVb', usecols = 2, dtype = 'float')
cont = np.loadtxt('NVb', usecols = 3, dtype = 'float')
cw= 1242.804
z_comp1 = 0.1385138493
cwz=(0.1384926549+1)*cw
cws=math.pow(cwz,2)
vel=[]
for i in range(len(lamda)):
    velocity = 300000*(((np.power(lamda[i],2))-cws)/((np.power(lamda[i],2))+cws))
    vel.append(velocity)
axs[0,1].plot(vel,flux,color='g',drawstyle='steps-mid', lw=0.7, alpha=0.8)
axs[0,1].plot(vel[50:82],cont[50:82],color='red',lw=1.2)
axs[0,1].plot(vel,error,color='y', drawstyle='steps-mid',lw=0.8)
axs[0,1].set_xlim(-1.900e+02,1.9e+02)
axs[0,1].set_ylim(-0.1e+00,1.3e+00)
axs[0,1].set_title('N V 1242',fontsize=5, x=0.2, y=0.2)
axs[0,1].tick_params(axis='both', which='major', labelsize=4.5)
axs[0,1].tick_params(axis='both', which='minor', labelsize=4.5)
axs[0,1].axhline(y=1.00,c="black",linewidth=0.6,ls='dotted')
axs[0,1].axhline(y=-0.1,c="black",linewidth=0.6,ls='dotted')
axs[0,1].axvline(0.5,c="black",linewidth=0.5,ls='dashed')
axs[0,1].spines['bottom'].set_visible(True)
axs[0,1].spines['top'].set_visible(True)
axs[0,1].text(65,0.4,'$λ_{obs}$ = 1414.949', fontsize=5)
axs[0,1].tick_params(left = False)



lamda = np.loadtxt('OVIa', usecols = 0, dtype = 'float')
flux  = np.loadtxt('OVIa', usecols = 1, dtype = 'float')
error = np.loadtxt('OVIa', usecols = 2, dtype = 'float')
cont = np.loadtxt('OVIa', usecols = 3, dtype = 'float')
cw= 1031.927
z_comp1 = 0.1385280701
cwz=(0.1384926549+1)*cw
cws=math.pow(cwz,2)
vel=[]
for i in range(len(lamda)):
    velocity = 300000*(((np.power(lamda[i],2))-cws)/((np.power(lamda[i],2))+cws))
    vel.append(velocity)
axs[1,1].plot(vel,flux,color='g',drawstyle='steps-mid', lw=0.7, alpha=0.8)
axs[1,1].plot(vel[50:100],cont[50:100],color='red',lw=1.2)
axs[1,1].plot(vel,error,color='y', drawstyle='steps-mid',lw=0.8)
axs[1,1].set_xlim(-1.900e+02,1.9e+02)
axs[1,1].set_ylim(-0.1e+00,1.3e+00)
axs[1,1].set_title('O VI 1032',fontsize=5, x=0.2, y=0.2)
axs[1,1].tick_params(axis='both', which='major', labelsize=4.5)
axs[1,1].tick_params(axis='both', which='minor', labelsize=4.5)
axs[1,1].axhline(y=1.00,c="black",linewidth=0.6,ls='dotted')
axs[1,1].axhline(y=-0.1,c="black",linewidth=0.6,ls='dotted')
axs[1,1].axvline(0.5,c="black",linewidth=0.5,ls='dashed')
axs[1,1].spines['bottom'].set_visible(True)
axs[1,1].spines['top'].set_visible(True)
axs[1,1].text(65,0.4,'$λ_{obs}$ =  1174.876', fontsize=5)
axs[1,1].tick_params(left = False)   


lamda = np.loadtxt('OVIb', usecols = 0, dtype = 'float')
flux  = np.loadtxt('OVIb', usecols = 1, dtype = 'float')
error = np.loadtxt('OVIb', usecols = 2, dtype = 'float')
cont = np.loadtxt('OVIb', usecols = 3, dtype = 'float')
cw= 1037.616
z_comp1 = 0.1385280701
cwz=(0.1384926549+1)*cw
cws=math.pow(cwz,2)
vel=[]
for i in range(len(lamda)):
    velocity = 300000*(((np.power(lamda[i],2))-cws)/((np.power(lamda[i],2))+cws))
    vel.append(velocity)
axs[2,1].plot(vel,flux,color='g',drawstyle='steps-mid', lw=0.7, alpha=0.8)
axs[2,1].plot(vel[50:80],cont[50:80],color='red',lw=1.2)
axs[2,1].plot(vel,error,color='y', drawstyle='steps-mid',lw=0.8)
axs[2,1].set_xlim(-1.900e+02,1.9e+02)
axs[2,1].set_ylim(-0.1e+00,1.3e+00)
axs[2,1].set_title('O VI 1038',fontsize=5, x=0.2, y=0.2)
axs[2,1].tick_params(axis='both', which='major', labelsize=4.5)
axs[2,1].tick_params(axis='both', which='minor', labelsize=4.5)
axs[2,1].axhline(y=1.00,c="black",linewidth=0.6,ls='dotted')
axs[2,1].axhline(y=-0.1,c="black",linewidth=0.6,ls='dotted')
axs[2,1].axvline(0.5,c="black",linewidth=0.5,ls='dashed')
axs[2,1].spines['bottom'].set_visible(True)
axs[2,1].spines['top'].set_visible(True)
axs[2,1].text(65,0.4,'$λ_{obs}$ = 1181.355', fontsize=5) 
axs[2,1].tick_params(left = False)

lamda = np.loadtxt('SiIIa', usecols = 0, dtype = 'float')
flux  = np.loadtxt('SiIIa', usecols = 1, dtype = 'float')
error = np.loadtxt('SiIIa', usecols = 2, dtype = 'float')
cont = np.loadtxt('SiIIa', usecols = 3, dtype = 'float')
cw= 1193.289
z_comp1 = 0.1385034484
cwz=(0.1384926549+1)*cw
cws=math.pow(cwz,2)
vel=[]
for i in range(len(lamda)):
    velocity = 300000*(((np.power(lamda[i],2))-cws)/((np.power(lamda[i],2))+cws))
    vel.append(velocity)
axs[3,1].plot(vel,flux,color='g',drawstyle='steps-mid', lw=0.7, alpha=0.8)
axs[3,1].plot(vel[50:76],cont[50:76],color='red',lw=1.2)
axs[3,1].plot(vel,error,color='y', drawstyle='steps-mid',lw=0.8)
axs[3,1].set_xlim(-1.900e+02,1.9e+02)
axs[3,1].set_ylim(-0.1e+00,1.3e+00)
axs[3,1].set_title('Si II 1193',fontsize=5, x=0.2, y=0.2)
axs[3,1].tick_params(axis='both', which='major', labelsize=4.5)
axs[3,1].tick_params(axis='both', which='minor', labelsize=4.5)
axs[3,1].axhline(y=1.00,c="black",linewidth=0.6,ls='dotted')
axs[3,1].axhline(y=-0.1,c="black",linewidth=0.6,ls='dotted')
axs[3,1].axvline(0.5,c="black",linewidth=0.5,ls='dashed')
axs[3,1].spines['bottom'].set_visible(True)
axs[3,1].spines['top'].set_visible(True)
axs[3,1].text(65,0.4,'$λ_{obs}$ = 1358.563', fontsize=5) 
axs[3,1].tick_params(left = False)

lamda = np.loadtxt('SiIIb', usecols = 0, dtype = 'float')
flux  = np.loadtxt('SiIIb', usecols = 1, dtype = 'float')
error = np.loadtxt('SiIIb', usecols = 2, dtype = 'float')
cont = np.loadtxt('SiIIb', usecols = 3, dtype = 'float')
cw= 1260.422
z_comp1 = 0.1385034484
cwz=(0.1384926549+1)*cw
cws=math.pow(cwz,2)
vel=[]
for i in range(len(lamda)):
    velocity = 300000*(((np.power(lamda[i],2))-cws)/((np.power(lamda[i],2))+cws))
    vel.append(velocity)
axs[4,1].plot(vel,flux,color='g',drawstyle='steps-mid', lw=0.7, alpha=0.8)
axs[4,1].plot(vel[50:87],cont[50:87],color='red',lw=1.2)
axs[4,1].plot(vel,error,color='y', drawstyle='steps-mid',lw=0.8)
axs[4,1].set_xlim(-1.900e+02,1.9e+02)
axs[4,1].set_ylim(-0.1e+00,1.3e+00)
axs[4,1].set_title('Si II 1260',fontsize=5, x=0.2, y=0.2)
axs[4,1].tick_params(axis='both', which='major', labelsize=4.5)
axs[4,1].tick_params(axis='both', which='minor', labelsize=4.5)
axs[4,1].axhline(y=1.00,c="black",linewidth=0.6,ls='dotted')
axs[4,1].axhline(y=-0.1,c="black",linewidth=0.6,ls='dotted')
axs[4,1].axvline(0.5,c="black",linewidth=0.5,ls='dashed')
axs[4,1].spines['bottom'].set_visible(True)
axs[4,1].spines['top'].set_visible(True)
axs[4,1].text(65,0.4,'$λ_{obs}$ = 1434.994', fontsize=5) 
axs[4,1].tick_params(left = False)

lamda = np.loadtxt('SiIII', usecols = 0, dtype = 'float')
flux  = np.loadtxt('SiIII', usecols = 1, dtype = 'float')
error = np.loadtxt('SiIII', usecols = 2, dtype = 'float')
cont = np.loadtxt('SiIII', usecols = 3, dtype = 'float')
cw= 1206.500
z_comp1 = 0.1384926549
cwz=(0.1384926549+1)*cw
cws=math.pow(cwz,2)
vel=[]
for i in range(len(lamda)):
    velocity = 300000*(((np.power(lamda[i],2))-cws)/((np.power(lamda[i],2))+cws))
    vel.append(velocity)
axs[5,1].plot(vel,flux,color='g',drawstyle='steps-mid', lw=0.7, alpha=0.8)
axs[5,1].plot(vel[50:99],cont[50:99],color='red',lw=1.2)
axs[5,1].plot(vel,error,color='y', drawstyle='steps-mid',lw=0.8)
axs[5,1].set_xlim(-1.900e+02,1.9e+02)
axs[5,1].set_ylim(-0.1e+00,1.3e+00)
axs[5,1].set_title('Si III 1206',fontsize=5, x=0.2, y=0.2)
axs[5,1].tick_params(axis='both', which='major', labelsize=4.5)
axs[5,1].tick_params(axis='both', which='minor', labelsize=4.5)
axs[5,1].axhline(y=1.00,c="black",linewidth=0.6,ls='dotted')
axs[5,1].axhline(y=-0.1,c="black",linewidth=0.6,ls='dotted')
axs[5,1].axvline(0.5,c="black",linewidth=0.5,ls='dashed')
axs[5,1].spines['bottom'].set_visible(True)
axs[5,1].spines['top'].set_visible(True)
axs[5,1].text(65,0.4,'$λ_{obs}$ = 1373.591', fontsize=5) 
axs[5,1].tick_params(left = False)

lamda = np.loadtxt('SiIVa', usecols = 0, dtype = 'float')
flux  = np.loadtxt('SiIVa', usecols = 1, dtype = 'float')
error = np.loadtxt('SiIVa', usecols = 2, dtype = 'float')
cont = np.loadtxt('SiIVa', usecols = 3, dtype = 'float')
cw= 1393.76018
z_comp1 = 0.138475707
cwz=(0.1384926549+1)*cw
cws=math.pow(cwz,2)
vel=[]
for i in range(len(lamda)):
    velocity = 300000*(((np.power(lamda[i],2))-cws)/((np.power(lamda[i],2))+cws))
    vel.append(velocity)
axs[6,1].plot(vel,flux,color='g',drawstyle='steps-mid', lw=0.7, alpha=0.8)
axs[6,1].plot(vel[51:75],cont[51:75],color='red',lw=1.2)
axs[6,1].plot(vel,error,color='y', drawstyle='steps-mid',lw=0.8)
axs[6,1].set_xlim(-1.900e+02,1.9e+02)
axs[6,1].set_ylim(-0.1e+00,1.3e+00)
axs[6,1].set_title('Si IV 1393',fontsize=5, x=0.2, y=0.2)
axs[6,1].tick_params(axis='both', which='major', labelsize=4.5)
axs[6,1].tick_params(axis='both', which='minor', labelsize=4.5)
axs[6,1].axhline(y=1.00,c="black",linewidth=0.6,ls='dotted')
axs[6,1].axhline(y=-0.1,c="black",linewidth=0.6,ls='dotted')
axs[6,1].axvline(0.5,c="black",linewidth=0.5,ls='dashed')
axs[6,1].spines['bottom'].set_visible(True)
axs[6,1].spines['top'].set_visible(True)
axs[6,1].text(65,0.4,'$λ_{obs}$ = 1586.762', fontsize=5)      
axs[6,1].tick_params(left = False)

lamda = np.loadtxt('SiIVb', usecols = 0, dtype = 'float')
flux  = np.loadtxt('SiIVb', usecols = 1, dtype = 'float')
error = np.loadtxt('SiIVb', usecols = 2, dtype = 'float')
cont = np.loadtxt('SiIVb', usecols = 3, dtype = 'float')
cw= 1402.77291
z_comp1 = 0.138475707
cwz=(0.1384926549+1)*cw
cws=math.pow(cwz,2)
vel=[]
for i in range(len(lamda)):
    velocity = 300000*(((np.power(lamda[i],2))-cws)/((np.power(lamda[i],2))+cws))
    vel.append(velocity)
axs[7,1].plot(vel,flux,color='g',drawstyle='steps-mid', lw=0.7, alpha=0.8)
axs[7,1].plot(vel[51:78],cont[51:78],color='red',lw=1.2)
axs[7,1].plot(vel,error,color='y', drawstyle='steps-mid',lw=0.8)
axs[7,1].set_xlim(-1.900e+02,1.9e+02)
axs[7,1].set_ylim(-0.1e+00,1.3e+00)
axs[7,1].set_title('Si IV 1402',fontsize=5, x=0.2, y=0.2)
axs[7,1].tick_params(axis='both', which='major', labelsize=4.5)
axs[7,1].tick_params(axis='both', which='minor', labelsize=4.5)
axs[7,1].axhline(y=1.00,c="black",linewidth=0.6,ls='dotted')
axs[7,1].axhline(y=-0.1,c="black",linewidth=0.6,ls='dotted')
axs[7,1].axvline(0.5,c="black",linewidth=0.5,ls='dashed')
axs[7,1].spines['bottom'].set_visible(True)
axs[7,1].spines['top'].set_visible(True)
axs[7,1].text(65,0.4,'$λ_{obs}$ = 1597.022', fontsize=5) 
axs[7,1].tick_params(left = False)


fig.supxlabel('Velocity ($km$ $s^{-1}$)', fontsize = 6)
fig.supylabel('Continuum Normalized Flux', fontsize = 6)
plt.subplots_adjust(left=0.08, bottom=0.08, right=0.98, top=0.99,wspace=0.00, hspace=0.0)
plt.show()
  
