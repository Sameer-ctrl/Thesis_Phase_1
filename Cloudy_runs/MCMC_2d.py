from astropy.io import fits
from astropy.table import Table
from numpy import *
from scipy.interpolate import interp2d
import emcee
import corner
import matplotlib.pyplot as plt
from roman import toRoman
from io import StringIO
import sys
import pickle

plt.style.use('../Python/my_style.mpl')

hdu=fits.open(f'Data/component_III_nH_Z_col_density_param.fits')
data=Table(hdu[1].data)

# ind=[]

# for i,row in enumerate(data):

#     if row['H']==0 or row['H']<10**16 or row['log_nH']>0 or row['log_Z']>0:
#        ind.append(i) 

# data.remove_rows([ind])


log_nH=data['log_nH']
log_Z=data['log_Z']

lognH_range=[-5,0]
logZ_range=[-3,2]

param=data['parameters']
ions=['Si+', 'Si+2','C+', 'C+2','O+5']

observations={'Si+':[13.19,0.41], 'Si+2':[12.87,0.08],'C+':[14.21,0.39], 'C+2':[13.81,0.04],'O+5':[13.91,0.04]}


def interp_func():

    func_dict={}
    # log_nH=data['log_nH']

    for i in observations.keys():

        x=[]

        for j in data[i].value:
            if j==0:
                x.append(10**(-30))
            
            else:
                x.append(j)

        log_col_den=log10(x)
        f=interp2d(log_nH,log_Z,log_col_den,kind='cubic')
        func_dict[i]=f

    return func_dict

interp_func_dict={}

for i in observations.keys():

    with open(f'Interp_2d_func/{i}_quintic.pkl','rb') as pickle_file:
        f=pickle.load(pickle_file)

    interp_func_dict[i]=f




# for i in arange(-5,0,0.5):

#     plt.scatter(arange(-3,0,0.01),interp_func_dict[ions[0]](i,arange(-3,0,0.01)),label=i)



# plt.legend()
# plt.show()
# quit()

# Si+: bad in middle and end
# Si+2: okayish, some waviness in middle bad in end
# C+: bad in end
# C+2: ok
# O+5: some ok some not


def log_posterior(theta,ions_to_use):

    lognH,logZ=theta

    #prior

    if lognH_range[0]<=lognH<=lognH_range[1] and logZ_range[0]<=logZ<=logZ_range[1]:
         log_prior=0
    
    else:
        log_prior=-inf

    #likelihood
    
    # interp_func_list=interp_func()
    model_col_den=[]
    observed_col_den=[]
    col_den_error=[]

    for i in ions_to_use:

        f=interp_func_dict[i]
        col_den=f(lognH,logZ)
        model_col_den.append(col_den)

        observed_col_den.append(observations[i][0])
        col_den_error.append(observations[i][1])
    
    model_col_den=array(model_col_den)    
    observed_col_den=array(observed_col_den)
    col_den_error=array(col_den_error)

    log_likelihood=-0.5*sum(log(2*pi*(col_den_error**2))+(((observed_col_den-model_col_den)/(sqrt(2)*col_den_error))**2))

    return log_prior+log_likelihood    # posterior


ions_to_use1=['Si+', 'Si+2','C+', 'C+2']
ions_to_use2=['Si+', 'Si+2','C+', 'C+2','O+5']

n_param = 2 
nwalkers = 50
nsteps = 5000

# n_guess = random.uniform(lognH_range[0],lognH_range[1], nwalkers)
# z_guess = random.uniform(logZ_range[0],logZ_range[1], nwalkers)

n_guess=random.normal(-1.6,0.2,nwalkers)
z_guess=random.normal(0.48,0.2,nwalkers)

starting_guesses = vstack((n_guess, z_guess)).T 

sampler1=emcee.EnsembleSampler(nwalkers, n_param, log_posterior, args=([ions_to_use1]))
sampler1.run_mcmc(starting_guesses, nsteps, progress=True)
# tau=sampler1.get_autocorr_time()  
# thin=int(mean(tau)/2) 
flat_samples1=sampler1.get_chain(flat=True)

n_guess=random.normal(-3.78,0.1,nwalkers)
z_guess=random.normal(-1.4,0.1,nwalkers)

# starting_guesses = vstack((n_guess, z_guess)).T 


sampler2=emcee.EnsembleSampler(nwalkers, n_param, log_posterior, args=([ions_to_use2]))
sampler2.run_mcmc(starting_guesses, nsteps, progress=True)
# tau=sampler2.get_autocorr_time()  
# thin=int(mean(tau)/2) 
flat_samples2=sampler2.get_chain(flat=True)

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

corner.corner(flat_samples1, labels=labels, quantiles=[0.16, 0.5, 0.84],show_titles=True, title_kwargs={"fontsize": 12},verbose=True,plot_contours=False)
quantiles1 = buffer.getvalue()
sys.stdout = sys.__stdout__
sol1=sol_write(quantiles1)


# fig.savefig('Files_n_figures/cornerplot_exc_OVI.png')

buffer = StringIO()
sys.stdout = buffer

corner.corner(flat_samples2, labels=labels, quantiles=[0.16, 0.5, 0.84],show_titles=True, title_kwargs={"fontsize": 12},verbose=True,plot_contours=False)
quantiles2 = buffer.getvalue()
sys.stdout = sys.__stdout__
sol2=sol_write(quantiles2)

# # fig.savefig('Files_n_figures/cornerplot_inc_OVI.png')
# plt.show()
# quit()

ions_all=[f'Si {toRoman(2)}', f'Si {toRoman(3)}',f'C {toRoman(2)}', f'C {toRoman(3)}',f'O {toRoman(6)}']
inds = random.randint(int(max(len(flat_samples1),len(flat_samples2))/1.33),min(len(flat_samples1),len(flat_samples2)), size=100)


plt.figure()

x=linspace(1,len(ions_all),len(ions_all))

for i in inds:

    sample1=flat_samples1[i]
    int_col_den1=[]
    int_col_den2=[]
    sample2=flat_samples2[i]

    for j in observations.keys():
        f=interp_func_dict[j]
        int_col_den1.append(f(sample1[0],sample1[1]))
        int_col_den2.append(f(sample2[0],sample2[1]))

    # plt.plot(x,int_col_den1,c='orange',alpha=0.1)
    if i!=inds[-1]:
        plt.plot(x,int_col_den1,c='orange',alpha=0.1)
        plt.plot(x,int_col_den2,c='green',alpha=0.1)
    
    else:
        plt.plot(x,int_col_den1,c='orange',alpha=0.4,label='Model samples (excluding O VI)')
        plt.plot(x,int_col_den2,c='green',alpha=0.4,label='Model samples (including O VI)')


# median_col_den_exc_OVI=array([interp_func_dict[i](sol1[0],sol1[3]) for i in observations.keys()])
# median_col_den_inc_OVI=array([interp_func_dict[i](sol2[0],sol2[3]) for i in observations.keys()])

median_col_den_exc_OVI=array([interp_func_dict[i](-1.6,0.48) for i in observations.keys()])
median_col_den_inc_OVI=array([interp_func_dict[i](-3.78,-1.40) for i in observations.keys()])

obs_col_den=array([observations[i][0]  for i in ions])
col_den_error=array([observations[i][1]  for i in ions])


chi_sq_exc=sum((((obs_col_den-median_col_den_exc_OVI)/col_den_error)**2)[:-1])
chi_sq_inc=sum((((obs_col_den-median_col_den_inc_OVI)/col_den_error)**2)[:-1])

print(f'reduced chi-sq excluding OVI : {chi_sq_exc/(4-2)}')
print(f'reduced chi-sq including OVI : {chi_sq_inc/(5-2)}')

plt.errorbar(x,obs_col_den,c='red',yerr=col_den_error, fmt='o',capsize=3,label='Observed')
plt.plot(x,median_col_den_exc_OVI,label='median solution (excluding O VI)',ls='--',lw=3,color='orange')
plt.plot(x,median_col_den_inc_OVI,label='median solution (including O VI)',ls='--',lw=3,color='green')
plt.xticks(x,ions_all,fontsize=20)
plt.yticks(fontsize=20)
plt.ylabel(r'$\mathbf{log \ (N \ {cm}^{-2})}$',labelpad=15)
plt.xlabel(r'$\mathbf{Ions}$',labelpad=15)
plt.legend()
# plt.savefig('Files_n_figures/Observed_and_predicted.png')
plt.show()





































# def obs_col_density(lognH,logZ,err=1):

#     for i in range(len(data)):
#         if log_nH[i]==lognH:
#             break
#         else:
#             pass

#     row=data[i]
#     obs_col_den=[]
#     obs_col_den_err=[]

#     for j in row:

#         j=log10(j)+logZ
#         sig=random.normal(0,err*0.01*j,1)
#         obs_col_den.append(j+sig)
#         obs_col_den_err.append(abs(sig))
    
#     return array(obs_col_den[:-2]),array(obs_col_den_err[:-2])





# def col_density(ion,*param_given):

#     "*param_given order: lognH, logZ, logT"

#     param_given=array(param_given)

#     for i in range(len(param)):

#         if array_equiv(param_given,param[i]):
#             break
#         else:
#             pass

#     return data[ion][i]

