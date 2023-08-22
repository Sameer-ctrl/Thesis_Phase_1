import numpy as np
from astropy.table import Table
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt
import emcee
import corner
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
# plt.rcParams.update({
#     "text.usetex": True,
#     "font.family": "ITC Bookman"
# })

#----data
def get_true_model(model_path, Q= 18, logZ = -1, uvb = 'KS19'):
    """
    :param model: The data where Q18 model is stored
    :return: a row of ion column densities at n_H = 1e-4 cm^-2
    """
    model = model_path + '/try_logZ({:.2f}).fits'.format(logZ)
    data = Table.read(model)
    true_ion_col = data [data['hden'] == 1e-4]
   # print(true_ion_col)
    return true_ion_col

#----model interpolation
def get_interp_func(model_path, ions_to_use, Q_uvb, uvb = 'KS19'):
    logZ = np.around(np.arange(-3, 0.05, 0.5), decimals = 2) # hardcoded
    #get nH array
    logZ_try = -1
    model_try = model_path + '/try_logZ({:.2f}).fits'.format(logZ_try)
    model = Table.read(model_try)
    lognH = np.log10(np.array(model['hden']))
    interpolation_function_list = []
    for ion in ions_to_use:
        z = np.zeros((len(lognH), len(logZ)))
        for i in range(len(logZ)):
            model = model_path + '/try_logZ({:.2f}).fits'.format(logZ[i])
            d = Table.read(model)
            d[ion][d[ion] == 0 ] = 1e-15 # for avoiding log10 (0) error
            z[:, i] = np.log10(d[ion]) #--- for log - log interpolation
        f = interp2d(lognH, logZ, z.T)
        interpolation_function_list.append(f)
    return interpolation_function_list

#----for mcmc
def log_likelihood(theta, interp_logf, obs_ion_col, col_err):
    """
    For a gaussian distributed errors
    :param theta: parameters [nH, Z]
    :param x: data x
    :param y: data y
    :param yerr: data err
    :return:
    """
    lognH, logZ=  theta
    # get metal ion column density for n_H and Z = 0.1
    col = []
    for i in range(len(obs_ion_col)):
        #print('==>', i, lognH, logT)
        #print(interp_logf[i](lognH, logT), i, lognH, logT)
        col_mod = interp_logf[i](lognH, logZ)[0]
        col.append(col_mod)
    model_col  = np.array(col)
    lnL = -0.5 * np.sum(np.log(2 * np.pi * col_err ** 2) + (obs_ion_col - model_col) ** 2 / col_err ** 2)
    return lnL

def log_prior(theta):
    lognH, logZ =  theta
    # flat prior
    if -5 <= lognH <= -2 and -3 <= logZ <= 0 :
        return 0.0
    return -np.inf

def log_posterior(theta, interp_func, data_col, sigma_col):
    log_p = log_prior(theta) + \
            log_likelihood(theta, interp_logf = interp_func, obs_ion_col = data_col, col_err = sigma_col)
    return log_p

def run_mcmc(model_path, data_col, sigma_col, Q_uvb, ions_to_use, uvb = 'KS19', figname = 'Corner_2D.pdf', same_error = False):
    # run_mcmc(model_Q= model, ions_to_use= ions)
    # ------------------ here is a way to run code
#    truths = [-4, -1]  # (lognH, logZ, logT) true values
    number_of_ions = len(ions_to_use)
#    data_col_all = get_true_model(model_path, Q=true_Q)
    # converting astropy table row to a list
#    data_col = []
#    for name in ions_to_use:
#        data_col.append(data_col_all[name][0])
#    np.random.seed(0)
#    if same_error:
#        sigma_col = 0.2 * np.ones(number_of_ions)
#    else:
#        sigma_col = np.random.uniform(0.01, 0.2, number_of_ions)
    print(data_col, sigma_col)
    interp_logf = get_interp_func(model_path = model_path, ions_to_use = ions_to_use, Q_uvb = Q_uvb, uvb = uvb)
    # Here we'll set up the computation. emcee combines multiple "walkers",
    # each of which is its own MCMC chain. The number of trace results will
    # be nwalkers * nsteps
    ndim = 2  # number of parameters in the model
    nwalkers = 50  # number of MCMC walkers
    nsteps = 10000  # number of MCMC steps to take
    # set theta near the maximum likelihood, with
    n_guess = np.random.uniform(-5, -2, nwalkers)
    z_guess = np.random.uniform(-3, 0, nwalkers)
    starting_guesses = np.vstack((n_guess, z_guess)).T  # initialise at a tiny sphere
    # Here's the function call where all the work happens:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=(interp_logf, data_col, sigma_col))
    sampler.run_mcmc(starting_guesses, nsteps, progress=True)
    # find out number of steps
    tau = sampler.get_autocorr_time()  # number of steps needed to forget the starting position
    #print(tau)
    thin = int(np.mean(tau) / 2)  # use this number for flattning the sample as done below
    #thin = 100
    full_flat_samples = sampler.get_chain(flat=True)
    flat_samples = sampler.get_chain(discard=thin * 20, flat=True)
    # we are discarding some initial steps roughly 5 times the autocorr_time steps
    # then we thin by about half the autocorrelation time steps for plotting => one does not have to do this step
    labels = [r'$\log [\textup{n}_\textup{H}/\textup{cm}^{-3}]$', r'$\log [\textup{Z}/\textup{Z}_\odot]$']
    #uvb_q= int((model_Q.split('try_Q')[-1]).split('.fits')[0])
    fig = corner.corner(flat_samples, labels=labels, fill_contours=True, quantiles=[0.16, 0.5, 0.84],show_titles=True, title_kwargs={"fontsize": 15},labelpad=-0.1)
    plt.savefig(figname,bbox_inches='tight')
    plt.close()
    for i in range(ndim):
        mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
        q = np.diff(mcmc)
        print(labels[i], '=', mcmc[1], q[0], q[1])
    return full_flat_samples,flat_samples, ndim, tau

ions_to_use= ['H','C+','Si+','N+2','O+5']
Q_uvb=18
model_path  = '/home/dheerajkumar/CLOUDY_code/0.39047/'
data_col=np.array([18.65,14.79,14.61,14.52,14.24])
sigma_col=np.array([0.20,0.05,0.08,0.06,0.05])

MCMC_1=run_mcmc(model_path, data_col, sigma_col, Q_uvb, ions_to_use, uvb = 'KS19', figname = '/home/dheerajkumar/CLOUDY_code/Corner_plot_1.pdf', same_error = False)