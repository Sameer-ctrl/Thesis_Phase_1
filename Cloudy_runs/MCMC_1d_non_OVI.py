from astropy.io import fits
from astropy.table import Table
from numpy import *
from scipy.interpolate import interp1d
import emcee
import corner
import matplotlib.pyplot as plt
from roman import toRoman
from io import StringIO
import sys

plt.style.use('../Python/my_style.mpl')


qso='he0056'
z_abs=0.043265
comp='II' 


ions=['C+3','Si+2','N+4']
col_den_dict=[[13.64,0.03],[12.37,0.07],[13.42,0.46]]

non_detc_ions=[]
non_detc_col_den=[]

N_Hi=15.30
logZ_ref=-1

observations=dict(zip(ions,col_den_dict))
non_detections=dict(zip(non_detc_ions,non_detc_col_den))

ions_to_use=array(ions)
# ions_to_use=array(ions)[[2,3,4,5]]


qso_list=loadtxt('../Python/Data/qso_list.txt',dtype=str)
qso_dict=dict(zip(qso_list[:,1],qso_list[:,0]))
qso_label=qso_dict[qso]

hdu=fits.open(f'{qso}/z={z_abs}/logZ={logZ_ref}/component_{comp}_PI_nH_col_density_param.fits')
data=Table(hdu[1].data)

log_nH=data['log_nH']
lognH_range=[-5,1]
logZ_range=[-4,2]


def ion_label(ion,ion_font_size=25,radicle_font_size=17):

    a=ion.split('+')

    if len(a)>1:

        if a[1]!='':
            return f'{{\\fontsize{{{ion_font_size}pt}}{{3em}}\selectfont{{}}$\mathbf{{{a[0]}}}$}} {{\\fontsize{{{radicle_font_size}pt}}{{3em}}\selectfont{{}}$\mathbf{{{toRoman(int(a[1])+1)}}}$}}'

        else:
            return f'{{\\fontsize{{{ion_font_size}pt}}{{3em}}\selectfont{{}}$\mathbf{{{a[0]}}}$}} {{\\fontsize{{{radicle_font_size}pt}}{{3em}}\selectfont{{}}$\mathbf{{{toRoman(2)}}}$}}'

    else:

        return f'{{\\fontsize{{{ion_font_size}pt}}{{3em}}\selectfont{{}}$\mathbf{{{a[0]}}}$}} {{\\fontsize{{{radicle_font_size}pt}}{{3em}}\selectfont{{}}$\mathbf{{{toRoman(1)}}}$}}'


def interp_func():

    func_dict={}

    for ion in list(observations.keys())+list(non_detections.keys()):

        ion_col_den=data[ion].value

        for i,j in enumerate(ion_col_den):
            if j==0:
                break

        log_nH_f=log_nH[:i]
        log_col_den=log10(ion_col_den[:i])
        f=interp1d(log_nH_f,log_col_den,fill_value='extrapolate',kind='quadratic')
        func_dict[ion]=f

    return func_dict

interp_func_dict=interp_func()


def log_posterior(theta,ions_to_use,non_detection_ion=non_detc_ions,prior_non_det=False):

    lognH,logZ=theta
        
    #prior
        
    if prior_non_det:  

        nd_bool=zeros(len(non_detection_ion))

        for j,n in enumerate(non_detection_ion):
            fn=interp_func_dict[n]
            col_den_nd=fn(lognH)+logZ-logZ_ref
            nd_bool[j]=col_den_nd<non_detections[n] 


        if lognH_range[0]<=lognH<=lognH_range[1] and logZ_range[0]<=logZ<=logZ_range[1] and all(nd_bool):
            log_prior=0

        else:
            log_prior=-inf

    else:

        if lognH_range[0]<=lognH<=lognH_range[1] and logZ_range[0]<=logZ<=logZ_range[1]:
            log_prior=0

        else:
            log_prior=-inf

    #likelihood
    
    model_col_den=[]
    observed_col_den=[]
    col_den_error=[]

    for i in ions_to_use:

        f=interp_func_dict[i]
        col_den=f(lognH)+logZ-logZ_ref
        model_col_den.append(col_den)

        observed_col_den.append(observations[i][0])
        col_den_error.append(observations[i][1])
    
    model_col_den=array(model_col_den)    
    observed_col_den=array(observed_col_den)
    col_den_error=array(col_den_error)

    log_likelihood=-0.5*sum(log(2*pi*(col_den_error**2))+(((observed_col_den-model_col_den)/(col_den_error))**2))

    return log_prior+log_likelihood    # posterior

def start_MCMC_samples(ions_to_use, guess=None, discard_tau=True, nwalkers=50, nsteps=5000,n_param = 2):

    if guess==None:

        n_guess = random.uniform(lognH_range[0],lognH_range[1], nwalkers)
        z_guess = random.uniform(logZ_range[0],logZ_range[1], nwalkers)
    
    else:

        if len(guess)==4:
            n_guess=random.normal(guess[0],guess[1],nwalkers)
            z_guess=random.normal(guess[2],guess[3],nwalkers)

        else:
            n_guess=ones(nwalkers)*float(guess[0])
            z_guess=ones(nwalkers)*float(guess[1])

    starting_guesses = vstack((n_guess, z_guess)).T 

    sampler=emcee.EnsembleSampler(nwalkers, n_param, log_posterior, args=([ions_to_use]))
    sampler.run_mcmc(starting_guesses, nsteps, progress=True)

    if discard_tau==True:
        tau=sampler.get_autocorr_time()  
        thin=int(mean(tau)/2) 
        flat_samples=sampler.get_chain(discard=thin*20,flat=True)
    
    else:
        flat_samples=sampler.get_chain(flat=True)

    return flat_samples


flat_samples=start_MCMC_samples(ions_to_use,discard_tau=True,nsteps=5000)


def sol_write(q):

    a=q.split(':')

    q1=a[1][2:-11]
    q2=a[2][2:-2]

    x=q1.split()
    y=q2.split()

    nH_quant=[float(x[1][:-2]),float(x[3][:-2]),float(x[5][:-2])]
    Z_quant=[float(y[1][:-2]),float(y[3][:-2]),float(y[5][:-2])]

    nH=[nH_quant[1],nH_quant[1]-nH_quant[0],nH_quant[2]-nH_quant[1]]
    Z=[Z_quant[1],Z_quant[1]-Z_quant[0],Z_quant[2]-Z_quant[1]]

    nH=[round(x,2) for x in nH]
    Z=[round(x,2) for x in Z]

    return nH+Z

labels=[r'$\mathbf{log \ n_H}$', r'$\mathbf{log \ Z}$']

buffer = StringIO()
sys.stdout = buffer

corner.corner(flat_samples, labels=labels, quantiles=[0.16, 0.5, 0.84],show_titles=True, title_kwargs={"fontsize": 12},verbose=True,plot_contours=False)
quantiles1 = buffer.getvalue()
sys.stdout = sys.__stdout__
sol=sol_write(quantiles1)

# fig.savefig('Files_n_figures/cornerplot_exc_OVI.png')


print(f'\nSolution : {sol}')


print('\n')
print(f'N(\ion{{H}}{{I}}) = {N_Hi}   \\\ \n')
print(f'Solution : $n_H$ = {sol[0]} $\pm$ {max([sol[1:3]])[0]} \hspace{{10mm}} $Z$ = {sol[3]} $\pm$ {max([sol[4:]])[0]} \\newline  \n')
# print(r'\newpage')
# print('')
# print(r'\textbf{Non-detections}')
# print('')
print(r'\begin{figure}[!htbp]')
print('    \centering')
print(f'    \includegraphics[width=1\linewidth]{{Ionisation-Modelling-Plots/{qso}-z={z_abs}-comp{comp}_logZ={logZ_ref}.png}}')
print(f'    \caption{{N(\ion{{H}}{{i}})={N_Hi}, log $Z_{{ref}}$={logZ_ref}}}')
print(r'\end{figure}')
# print(f'\includegraphics[width=1\linewidth]{{Ionisation-Modelling-Plots/{qso}-z={z_abs}-comp{comp}_logZ={logZ_ref}.png}}\n')

inds = random.randint(int(len(flat_samples)/1.33),len(flat_samples), size=100)

plt.figure(figsize=(16,9))

ions_roman=[ion_label(i) for i in list(observations.keys())+list(non_detections.keys())]
x=linspace(1,len(ions_roman),len(ions_roman))

for i in inds:

    sample=flat_samples[i]
    int_col_den=[]

    for j in list(observations.keys())+list(non_detections.keys()):
        f=interp_func_dict[j]
        int_col_den.append(f(sample[0])+sample[1]-logZ_ref)

    if i!=inds[-1]:
        plt.plot(x,int_col_den,c='orange',alpha=0.1)
    
    else:
        plt.plot(x,int_col_den,c='orange',alpha=0.4,label=f'Model samples')

median_col_den=array([interp_func_dict[i](sol[0])+sol[3]-logZ_ref for i in list(observations.keys())+list(non_detections.keys())])

obs_col_den=array([observations[i][0]  for i in ions])
col_den_error=array([observations[i][1]  for i in ions])


chi_sq=sum((((obs_col_den-median_col_den[:len(ions)])/col_den_error)**2)[:-1])

print(f'\nreduced chi-sq: {chi_sq/(4-2)}')


plt.errorbar(x[:len(ions)],obs_col_den,c='red',yerr=col_den_error, fmt='o',capsize=3,label='Observed')
# plt.errorbar(x[len(ions):],non_detc_col_den,yerr=0.25,xerr=0.1,uplims=1,fmt='o',label='Non-detections')
plt.plot(x,median_col_den,label=f'median solution',ls='--',lw=3,color='orange')
plt.xticks(x,ions_roman,fontsize=20)
plt.yticks(fontsize=20)
plt.ylabel(r'$\mathbf{log \ (N \ {cm}^{-2})}$',labelpad=15)
plt.xlabel(r'$\mathbf{Ions}$',labelpad=15)
plt.legend()
plt.title(f'$\mathbf{{{qso_label} \ (z_{{abs}}={z_abs})}}$',fontsize=30)
plt.savefig(f'../LaTeX/BLA_Survey_results/Ionisation-Modelling-Plots/{qso}-z={z_abs}-comp{comp}_logZ={logZ_ref}.png')
plt.show()
