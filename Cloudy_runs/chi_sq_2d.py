from operator import indexOf
from astropy.io import fits,ascii
from astropy.table import Table
from numpy import *
from roman import toRoman
from scipy.interpolate import interp2d,interp1d
import pickle
import matplotlib.pyplot as plt

import warnings

warnings.filterwarnings("ignore")

plt.style.use('Files_n_figures/my_style.mpl')

# def ion_label(ion,state):

#     # return f'{{\\fontsize{{25pt}}{{3em}}\selectfont{{}}$\mathbf{{{ion}}}$}} {{\\fontsize{{17pt}}{{3em}}\selectfont{{}}$\mathbf{{{state}}}$}}'

def ion_label(ion,ion_font_size=25,radicle_font_size=17):

    a=ion.split('+')

    if len(a)>1:

        if a[1]!='':
            return f'{{\\fontsize{{{ion_font_size}pt}}{{3em}}\selectfont{{}}$\mathbf{{{a[0]}}}$}} {{\\fontsize{{{radicle_font_size}pt}}{{3em}}\selectfont{{}}$\mathbf{{{toRoman(int(a[1])+1)}}}$}}'

        else:
            return f'{{\\fontsize{{{ion_font_size}pt}}{{3em}}\selectfont{{}}$\mathbf{{{a[0]}}}$}} {{\\fontsize{{{radicle_font_size}pt}}{{3em}}\selectfont{{}}$\mathbf{{{toRoman(2)}}}$}}'

    else:

        return f'{{\\fontsize{{{ion_font_size}pt}}{{3em}}\selectfont{{}}$\mathbf{{{a[0]}}}$}} {{\\fontsize{{{radicle_font_size}pt}}{{3em}}\selectfont{{}}$\mathbf{{{toRoman(1)}}}$}}'





ions=['Si+2', 'Si+','C+2', 'C+','O+5']
ions_roman=[ion_label('Si+2',ion_font_size=35,radicle_font_size=23),ion_label('Si+',ion_font_size=35,radicle_font_size=23),ion_label('C+2',ion_font_size=35,radicle_font_size=23),ion_label('C+',ion_font_size=35,radicle_font_size=23),ion_label('O+5',ion_font_size=35,radicle_font_size=23)]


# observations={'Si+2':[12.87,0.08],'Si+':[13.19,0.41], 'C+2':[13.72,0.04],'C+':[14.21,0.39], 'O+5':[13.91,0.04]}
observations={'Si+2':[12.90,0.02+0.15],'Si+':[12.90,0.16+0.18], 'C+2':[13.88,0.01+0.09],'C+':[13.65,0.19+0.14], 'O+5':[13.92,0.01+0.05]}


interp_func_dict={}

for i in observations.keys():

    with open(f'Interp_2d_func_Danforth/{i}_quintic.pkl','rb') as pickle_file:
        f=pickle.load(pickle_file)
    
    interp_func_dict[i]=f

obs_col_den=array([observations[i][0]  for i in ions])
col_den_error=array([observations[i][1]  for i in ions]) 

def chi_sq_func(nH,Z):

    if type(nH)==type(1) or type(nH)==type(1.0) or type(nH)==type(float64(1)):
        nH=[nH]

    else:
        nH=nH

    col_den=array([interp_func_dict[i](nH,Z) for i in observations.keys()])
    chi_sq_exc=zeros(len(nH))
    chi_sq_inc=zeros(len(nH))

    for i in range(len(nH)):

        med_col_den=col_den[:,i]
        chi_sq_exc[i]=sum((((obs_col_den-med_col_den)/col_den_error)**2)[:-1])
        chi_sq_inc[i]=sum((((obs_col_den-med_col_den)/col_den_error)**2))

    return chi_sq_exc,chi_sq_inc


nH_global=arange(-5,0.01,0.01)
Z_global=arange(-3,2.01,0.01)

x=array(list(nH_global)*len(Z_global))
y=zeros(len(x))

chi_sq_exc=zeros((len(Z_global),len(nH_global)))
chi_sq_inc=zeros((len(Z_global),len(nH_global)))

for i,a in enumerate(Z_global):

    chi_sq_exc[i],chi_sq_inc[i]=chi_sq_func(nH_global,a)
    y[len(nH_global)*i:len(nH_global)*(i+1)]=a
observations={'Si+2':[12.87,0.08],'Si+':[13.19,0.41], 'C+2':[13.72,0.04],'C+':[14.21,0.39], 'O+5':[13.91,0.04]}


def model(chi_sq,name,c):

    chi_sq_arr=chi_sq.flatten()
    ind=argmin(chi_sq_arr)

    nH=round(x[ind],2)
    Z=round(y[ind],2)
    chi_sq_min=round(chi_sq_arr[ind],3)

    if name[:3]=='Exc':
        a=0
    else:
        a=1
    
    chi_sq_err_nH=[chi_sq_func(n,Z)[a] for n in nH_global]
    chi_sq_err_Z=[chi_sq_func(nH,z)[a] for z in Z_global]

    nH_level=[]
    Z_level=[]

    i=indexOf(nH_global,x[ind])

    for j in range(i):
        if chi_sq_err_nH[i-j]-chi_sq_min>1:
            nH_level.append(round(nH_global[i-j],2))
            break

    for j in range(i,len(nH_global)):
        if chi_sq_err_nH[j]-chi_sq_min>1:
            nH_level.append(round(nH_global[j],2))
            break

    i=indexOf(Z_global,y[ind])


    for j in range(i):
        if chi_sq_err_Z[i-j]-chi_sq_min>1:
            Z_level.append(round(Z_global[i-j],2))
            break

    for j in range(i,len(Z_global)):
        if chi_sq_err_Z[j]-chi_sq_min>1:
            Z_level.append(round(Z_global[j],2))
            break

    nH_err=[round(nH-nH_level[0],2),round(nH_level[1]-nH,2)]
    Z_err=[round(Z-Z_level[0],2),round(Z_level[1]-Z,2)]

    print(f'{name} : nH = {nH}  [{nH_err[0]},{nH_err[1]}]   Z = {Z}  [{Z_err[0]},{Z_err[1]}]  chi-sq = {chi_sq_min}')

    col_den=array([interp_func_dict[ion](nH,Z) for ion in observations.keys()])
    xaxis=linspace(1,len(ions_roman),len(ions_roman))

    plt.plot(xaxis,col_den,ls='--',lw=3,color=c)
    plt.scatter(xaxis,col_den,color=c,s=100,label=name)

    return float(nH),float(Z),nH_err,Z_err


def plot_samples(m,c,n=50,i1=1,i2=1):

    nH_sample=random.normal(m[0],max(m[2])*i1,n)
    Z_sample=random.normal(m[1],max(m[3])*i2,n)

    sample = vstack((nH_sample, Z_sample)).T 
    xaxis=linspace(1,len(ions_roman),len(ions_roman))

    for s in sample:
        nH,Z=s
        col_den=array([interp_func_dict[i](round(nH,3),round(Z,3)) for i in observations.keys()]) 
        plt.plot(xaxis,col_den,alpha=0.05,color=c)


xaxis=linspace(1,len(ions_roman),len(ions_roman))

plt.figure(figsize=(16,10),dpi=300)

m1=model(chi_sq_exc,'Excluding OVI','orange')
m2=model(chi_sq_inc,'Including OVI','green')
plt.clf()

plt.errorbar(xaxis,obs_col_den,c='red',yerr=col_den_error, fmt='o',capsize=3,label=r'$\mathbf{Observed}$')
plot_samples(m1,'orange',n=100)
m1=model(chi_sq_exc,r'$\mathbf{Excluding \ }$'+ion_label('O+5'),'orange')
plot_samples(m2,'green',n=100)
m2=model(chi_sq_inc,r'$\mathbf{Including \ }$'+ion_label('O+5'),'green')
# plt.text(1,9,r'Excluding OVI : log nH = -2.24$\pm$0.03 \ \ \ \ \ \  log Z = -0.31$\pm$0.06 \ \ \ \ \ \ $\chi^{2}=4.268$')
# plt.text(1,8.5,r'Including OVI : log nH = -3.88$\pm$0.02 \ \ \ \ \ \ log Z = -1.51$\pm$0.03 \ \ \ \ \ \ $\chi^{2}=275.666$')
plt.xticks(xaxis,ions_roman,fontsize=30)
plt.yticks([8,9,10,11,12,13,14,15],fontsize=30)
plt.ylabel(r'$\mathbf{log \ [N \ ({cm}^{-2})]}$',labelpad=25,fontsize=40)
plt.xlabel(r'$\mathbf{Ions}$',labelpad=18,fontsize=40)
# plt.title(r'Solution using $\chi^{2}$ minimization',pad=15)
plt.legend(loc='lower left',fontsize=25)
plt.ylim(7.11,15.28)

plt.savefig('Observed_and_predicted_Danforth.png')

# plt.show()





# fig=plt.figure()
# ax=plt.axes(projection ='3d')

# ax.scatter(x,y,z)
# ax.set_xlabel('nH')
# ax.set_ylabel('Z')
# ax.set_zlabel('chi_square')

























