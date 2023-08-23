import numpy as np
import astropy.table as tab
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt
import emcee
import corner

def get_interp_func(model_path, ions_to_use, Q_uvb, uvb = 'KS19'):
    logZ = np.around(np.arange(-3, 0.05, 1), decimals = 2) # hardcoded
    #get nH array
    logZ_try = -1
    model_try = model_path + '/try_logZ({:.2f}).fits'.format(logZ_try)
    model = tab.Table.read(model_try)
    lognH = np.log10(np.array(model['hden']))
    interpolation_function_list = []
    for ion in ions_to_use:
        z = np.zeros((len(lognH), len(logZ)))
        for i in range(len(logZ)):
            model = model_path + '/try_logZ({:.2f}).fits'.format(logZ[i])
            d = tab.Table.read(model)
            d[ion][d[ion] == 0 ] = 1e-15 # for avoiding log10 (0) error
            z [:, i] = np.log10(d[ion]) #--- for log - log interpolation
        f = interp1d(lognH, z.T, fill_value='extrapolate')
        interpolation_function_list.append(f)
    return interpolation_function_list

def log_likelihood(theta, interp_logf, obs_ion_col, col_err, reference_log_metal = -2.0):
    """
    For a gaussian distributed errors
    :param theta: parameters [nH, Z]
    :param x: data x
    :param y: data y
    :param yerr: data err
    :return:
    """
    lognH, logZ =  theta
    # get metal ion column density for n_H and Z = 0.1
    col = []
    for i in range(len(obs_ion_col)):
        #print('==>', i, lognH, logT)
        #print(interp_logf[i](lognH, logT), i, lognH, logT)
        col_mod = 10**interp_logf[i](lognH)[1]
        col.append(col_mod)
    # scale the column densities by the metallicity Z
    metal_scaling_linear = 10 ** logZ / 10 ** reference_log_metal
    model_col = np.log10(np.clip(col, 1e-10, 1e22) * metal_scaling_linear)
    lnL = -0.5 * np.sum(np.log(2 * np.pi * col_err ** 2) + (obs_ion_col - model_col) ** 2 / col_err ** 2)
    return lnL

def log_prior(theta):
    lognH, logZ=  theta
    # flat prior
    if -5.0 < lognH < -2.0 and -3.0 < logZ < 0.0:   # Better result when ranges are same as grids.
        return 0.0
    return -np.inf

def log_posterior(theta, interp_func, data_col, sigma_col):
    log_p = log_prior(theta) + \
            log_likelihood(theta, interp_logf = interp_func, obs_ion_col = data_col, col_err = sigma_col)
    return log_p

def run_mcmc(model_path, data_col, sigma_col, Q_uvb, ions_to_use, uvb = 'KS19', figname = 'Corner.pdf'):
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
    labels = ['log nH', 'log Z']
    #uvb_q= int((model_Q.split('try_Q')[-1]).split('.fits')[0])
    fig = corner.corner(flat_samples, labels=labels, quantiles=[0.16, 0.5, 0.84],show_titles=True, title_kwargs={"fontsize": 12})
    fig.savefig(figname)
    plt.close()
    for i in range(ndim):
        mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
        q = np.diff(mcmc)
        print(labels[i], '=', mcmc[1], q[0], q[1])
    return full_flat_samples,flat_samples, ndim, tau

ions_to_use= ['H','C+','Si+','N+2','O+5']
Q_uvb=18
model_path  = '/home/dheerajkumar/CLOUDY_code/'
data_col=np.array([18.65,16.635274,14.650345,14.729093,14.240797])
sigma_col=np.array([0.2,1.414981,0.058193,0.077588,0.024540])
MCMC=run_mcmc(model_path, data_col, sigma_col, Q_uvb, ions_to_use, uvb = 'KS19', figname = '/home/dheerajkumar/CLOUDY_code/Corner_1D_1.pdf')