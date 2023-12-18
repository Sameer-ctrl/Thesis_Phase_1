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


qso='s135712'
comp='III'
z_abs=0.097869

ions=['Si+2','Si+3','C+','C+3','O+5']
col_den_dict=[[14.67,0.13],[13.3,0.05],[13.92,0.04],[13.63,0.12],[14.3,0.11]]

observations=dict(zip(ions,col_den_dict))
# observations={'C+3':[13.58,0.09],'Si+2':[12.58,0.05],'Si+3':[12.69,0.1],'O+5':[13.77,0.11]}

ions_to_use1=ions[:-1]
ions_to_use2=ions

qso_list=loadtxt('../Python/Data/qso_list.txt',dtype=str)
qso_dict=dict(zip(qso_list[:,1],qso_list[:,0]))
qso_label=qso_dict[qso]

hdu=fits.open(f'{qso}/z={z_abs}/component_{comp}_PI_nH_col_density_param.fits')
data=Table(hdu[1].data)

log_nH=data['log_nH']
logZ_ref=-1
lognH_range=[-5,1]
logZ_range=[-3,2]


def ion_label(ion,ion_font_size=25,radicle_font_size=17):

    a=ion.split('+')

    if a[1]!='':
        return f'{{\\fontsize{{{ion_font_size}pt}}{{3em}}\selectfont{{}}$\mathbf{{{a[0]}}}$}} {{\\fontsize{{{radicle_font_size}pt}}{{3em}}\selectfont{{}}$\mathbf{{{toRoman(int(a[1])+1)}}}$}}'

    else:
        return f'{{\\fontsize{{{ion_font_size}pt}}{{3em}}\selectfont{{}}$\mathbf{{{a[0]}}}$}} {{\\fontsize{{{radicle_font_size}pt}}{{3em}}\selectfont{{}}$\mathbf{{{toRoman(2)}}}$}}'


def interp_func():

    func_dict={}

    for ion in observations.keys():

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


def log_posterior(theta,ions_to_use):

    lognH,logZ=theta

    #prior

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


flat_samples_exc_OVI=start_MCMC_samples(ions_to_use1,discard_tau=True,nsteps=5000)
flat_samples_inc_OVI=start_MCMC_samples(ions_to_use2,discard_tau=True,nsteps=5000)


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

labels=['log nH', 'log Z']

buffer = StringIO()
sys.stdout = buffer

corner.corner(flat_samples_exc_OVI, labels=labels, quantiles=[0.16, 0.5, 0.84],show_titles=True, title_kwargs={"fontsize": 12},verbose=True,plot_contours=False)
quantiles1 = buffer.getvalue()
sys.stdout = sys.__stdout__
sol1=sol_write(quantiles1)

# fig.savefig('Files_n_figures/cornerplot_exc_OVI.png')

buffer = StringIO()
sys.stdout = buffer

corner.corner(flat_samples_inc_OVI, labels=labels, quantiles=[0.16, 0.5, 0.84],show_titles=True, title_kwargs={"fontsize": 12},verbose=True,plot_contours=False)
quantiles2 = buffer.getvalue()
sys.stdout = sys.__stdout__
sol2=sol_write(quantiles2)

# fig.savefig('Files_n_figures/cornerplot_inc_OVI.png')

print(f'\nExcluding OVI : {sol1}')
print(f'Including OVI : {sol2}\n')

print(f'Excluding \ion{{O}}{{vi}} : $n_H$ = {sol1[0]} $\pm$ {max([sol1[1:3]])[0]} \hspace{{10mm}} $Z$ = {sol1[3]} $\pm$ {max([sol1[4:]])[0]}\n')
print(f'Including \ion{{O}}{{vi}} : $n_H$ = {sol2[0]} $\pm$ {max([sol2[1:3]])[0]} \hspace{{10mm}} $Z$ = {sol2[3]} $\pm$ {max([sol2[4:]])[0]}')
print(f'\\\\\\\\')
print(f'\includegraphics[width=1\linewidth]{{Ionisation-Modelling-Plots/{qso}-z={z_abs}-comp{comp}.png}}\n')

inds = random.randint(int(max(len(flat_samples_exc_OVI),len(flat_samples_inc_OVI))/1.33),min(len(flat_samples_exc_OVI),len(flat_samples_inc_OVI)), size=100)

plt.figure(figsize=(16,9))

ions_roman=[ion_label(i) for i in ions]
x=linspace(1,len(ions_roman),len(ions_roman))

for i in inds:

    sample1=flat_samples_exc_OVI[i]
    int_col_den1=[]
    int_col_den2=[]
    sample2=flat_samples_inc_OVI[i]

    for j in observations.keys():
        f=interp_func_dict[j]
        int_col_den1.append(f(sample1[0])+sample1[1]-logZ_ref)
        int_col_den2.append(f(sample2[0])+sample2[1]-logZ_ref)

    if i!=inds[-1]:
        plt.plot(x,int_col_den1,c='orange',alpha=0.1)
        plt.plot(x,int_col_den2,c='green',alpha=0.1)
    
    else:
        plt.plot(x,int_col_den1,c='orange',alpha=0.4,label=f'Model samples (excluding {ion_label("O+5",ion_font_size=15,radicle_font_size=10)})')
        plt.plot(x,int_col_den2,c='green',alpha=0.4,label=f'Model samples (including {ion_label("O+5",ion_font_size=15,radicle_font_size=10)})')


median_col_den_exc_OVI=array([interp_func_dict[i](sol1[0])+sol1[3]-logZ_ref for i in observations.keys()])
median_col_den_inc_OVI=array([interp_func_dict[i](sol2[0])+sol2[3]-logZ_ref for i in observations.keys()])

obs_col_den=array([observations[i][0]  for i in ions])
col_den_error=array([observations[i][1]  for i in ions])


chi_sq_exc=sum((((obs_col_den-median_col_den_exc_OVI)/col_den_error)**2)[:-1])
chi_sq_inc=sum((((obs_col_den-median_col_den_inc_OVI)/col_den_error)**2)[:-1])

print(f'\nreduced chi-sq excluding OVI : {chi_sq_exc/(4-2)}')
print(f'reduced chi-sq including OVI : {chi_sq_inc/(5-2)}')


plt.errorbar(x,obs_col_den,c='red',yerr=col_den_error, fmt='o',capsize=3,label='Observed')
plt.plot(x,median_col_den_exc_OVI,label=f'median solution (excluding {ion_label("O+5",ion_font_size=15,radicle_font_size=10)})',ls='--',lw=3,color='orange')
plt.plot(x,median_col_den_inc_OVI,label=f'median solution (including {ion_label("O+5",ion_font_size=15,radicle_font_size=10)})',ls='--',lw=3,color='green')
plt.xticks(x,ions_roman,fontsize=20)
plt.yticks(fontsize=20)
plt.ylabel(r'$\mathbf{log \ (N \ {cm}^{-2})}$',labelpad=15)
plt.xlabel(r'$\mathbf{Ions}$',labelpad=15)
plt.legend()
plt.title(f'$\mathbf{{{qso_label} \ (z_{{abs}}={z_abs})}}$',fontsize=30)
plt.savefig(f'../LaTeX/BLA Survey results/Ionisation-Modelling-Plots/{qso}-z={z_abs}-comp{comp}.png')
plt.show()
