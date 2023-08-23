from astropy.io import fits
from astropy.table import Table
from numpy import *
from scipy.interpolate import interp1d
import emcee
import corner
import matplotlib.pyplot as plt

# plt.style.use('../Python/my_style.mpl')

hdu=fits.open('Data/component_III_1d_col_density_param.fits')
data=Table(hdu[1].data)

mask=data['H']!=0
data=data[mask]

log_nH=data['nH']
logZ_ref=-1

lognH_range=[-6,2]
logZ_range=[-3,2]

param=data['parameters']
ions=['C+2', 'C+3', 'O+6', 'Si+2', 'Si+3']


def interp_func():

    func_list=[]
    log_nH=data['nH']

    for i in ions:
        log_col_den=log10([data[i].value])
        f=interp1d(log_nH,log_col_den,fill_value='extrapolate')
        func_list.append(f)

    return func_list


def obs_col_density(lognH,logZ,err=1):

    for i in range(len(data)):
        if log_nH[i]==lognH:
            break
        else:
            pass

    row=data[i]
    obs_col_den=[]
    obs_col_den_err=[]

    for j in row:
        j=log10(j)+logZ-logZ_ref
        sig=random.normal(0,err*0.01*j,1)
        obs_col_den.append(j+sig)
        obs_col_den_err.append(abs(sig))
    
    return array(obs_col_den[:-2]),array(obs_col_den_err[:-2])


def log_posterior(theta,observed_col_den,col_den_err):

    lognH,logZ=theta

    #prior

    if lognH_range[0]<=lognH<=lognH_range[1] and logZ_range[0]<=logZ<=logZ_range[1]:
         log_prior=0
    
    else:
        log_prior=-inf

    #likelihood
    
    interp_func_list=interp_func()
    model_col_den=[]

    for f in interp_func_list:
        col_den=f(lognH)+logZ-logZ_ref
        model_col_den.append(col_den)
    
    model_col_den=array(model_col_den)

    log_likelihood=-0.5*sum(log(2*pi*(col_den_err**2))+(((observed_col_den-model_col_den)/(sqrt(2)*col_den_err))**2))

    return log_prior+log_likelihood    # posterior


# obs=obs_col_density(-4.68,-1.25,err=4)
# ions=['C+2', 'C+3', 'O+6', 'Si+2', 'Si+3']

obs_col_density=[13.49,13.65,13.84,12.53,12.69]
col_density_error=[0.05,0.03,0.03,0.07,0.06]

obs=[obs_col_density,col_density_error]


print(f'obs column density : {obs[0].T[0]} \n')
print(f'column density err : {obs[1].T[0]} \n')
print(f'% error : {transpose((obs[1]/obs[0])*100)[0]} \n')

n_param = 2 
nwalkers = 50
nsteps = 10000  
n_guess = random.uniform(lognH_range[0],lognH_range[1], nwalkers)
z_guess = random.uniform(logZ_range[0],logZ_range[1], nwalkers)
starting_guesses = vstack((n_guess, z_guess)).T 

sampler=emcee.EnsembleSampler(nwalkers, n_param, log_posterior, args=(obs[0],obs[1]))
sampler.run_mcmc(starting_guesses, nsteps, progress=True)

tau=sampler.get_autocorr_time()  
thin=int(mean(tau)/2) 

flat_samples=sampler.get_chain(discard=thin * 20, flat=True)
labels=['log nH', 'log Z']

fig=corner.corner(flat_samples, labels=labels, quantiles=[0.16, 0.5, 0.84],show_titles=True, title_kwargs={"fontsize": 12})
fig.savefig('Files_n_figures/cornerplot.png')

plt.show()































# def col_density(ion,*param_given):

#     "*param_given order: lognH, logZ, logT"

#     param_given=array(param_given)

#     for i in range(len(param)):

#         if array_equiv(param_given,param[i]):
#             break
#         else:
#             pass

#     return data[ion][i]

