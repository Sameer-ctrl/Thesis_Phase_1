from astropy.io import fits
from astropy.table import Table
from matplotlib import markers
from numpy import *
from scipy.interpolate import interp1d
import emcee
import corner
import matplotlib.pyplot as plt
from roman import toRoman
from io import StringIO
import sys

plt.style.use('../Python/my_style.mpl')

# observations={'Si+2':[13.19,0.21], 'Si+3':[12.87,0.19],'C+2':[14.21,0.21], 'C+3':[13.81,0.23],'O+6':[13.91,0.22]}

#rebinned observation
observations={'Si+':[13.19,0.41], 'Si+2':[12.87,0.08],'C+':[14.21,0.39], 'C+2':[13.81,0.04],'O+5':[13.91,0.04]}

n=7

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


def run_UVB(ions_to_use,file_name):

    sol=[]
    interp_func_list=[]

    for q in range(14,14+n):

        print(f'UV_q = {q}')

        hdu=fits.open(f'Data/UV_Q_new/component_III_nH_UV_Q{q}_col_density_param.fits')
        data=Table(hdu[1].data)
        print

        mask=data['H']!=0
        data=data[mask]

        log_nH=data['log_nH']

        logZ_ref=-1

        lognH_range=[-5,2]
        logZ_range=[-3,2]

        # param=data['parameters']
        # ions=['Si+2', 'Si+3','C+2', 'C+3','O+6']

        # observations={'Si+2':[12.53,0.07], 'Si+3':[12.69,0.06],'C+2':[13.49,0.05], 'C+3':[13.65,0.03],'O+6':[13.84,0.03]}
        
        #rebinned observation
        

        def interp_func():

            func_dict={}
            # log_nH=data['log_nH']

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

        # for i in ions:
        #     plt.plot(log_nH,interp_func_dict[i](log_nH),label=i)

        # plt.legend()
        # plt.show()


        # quit()

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
                col_den=f(lognH)+logZ-logZ_ref
                model_col_den.append(col_den)

                observed_col_den.append(observations[i][0])
                col_den_error.append(observations[i][1])
            
            model_col_den=array(model_col_den)    
            observed_col_den=array(observed_col_den)
            col_den_error=array(col_den_error)

            log_likelihood=-0.5*sum(log(2*pi*(col_den_error**2))+(((observed_col_den-model_col_den)/(sqrt(2)*col_den_error))**2))

            return log_prior+log_likelihood    # posterior



        # ions_to_use=['Si+2', 'Si+3','C+2', 'C+3']
        # ions_to_use2=['Si+2', 'Si+3','C+2', 'C+3','O+6']

        n_param = 2 
        nwalkers = 50
        nsteps = 3000  
        # n_guess = random.uniform(lognH_range[0],lognH_range[1], nwalkers)
        # z_guess = random.uniform(logZ_range[0],logZ_range[1], nwalkers)

        n_guess=random.normal(-3.8,0.3,nwalkers)
        z_guess=random.normal(-0.9,0.2,nwalkers)

        starting_guesses = vstack((n_guess, z_guess)).T 

        sampler1=emcee.EnsembleSampler(nwalkers, n_param, log_posterior, args=([ions_to_use]))
        sampler1.run_mcmc(starting_guesses, nsteps, progress=True)
        tau=sampler1.get_autocorr_time()  
        thin=int(mean(tau)/2) 
        flat_samples1=sampler1.get_chain(discard=thin*20, flat=True)

        # n_guess=random.normal(-3.8,0.2,nwalkers)
        # z_guess=random.normal(-1.5,0.2,nwalkers)

        # starting_guesses = vstack((n_guess, z_guess)).T 


        # sampler2=emcee.EnsembleSampler(nwalkers, n_param, log_posterior, args=([ions_to_use2]))
        # sampler2.run_mcmc(starting_guesses, nsteps, progress=True)
        # tau=sampler2.get_autocorr_time()  
        # thin=int(mean(tau)/2) 
        # flat_samples2=sampler2.get_chain(discard=thin*20, flat=True)


        labels=['log nH', 'log Z']
        # print('Excluding O VI \n')
        buffer = StringIO()
        sys.stdout = buffer
        fig1=corner.corner(flat_samples1, labels=labels, quantiles=[0.16, 0.5, 0.84],show_titles=True, title_kwargs={"fontsize": 12},verbose=True)
        quantiles1 = buffer.getvalue()
        sys.stdout = sys.__stdout__
        med_sol=sol_write(quantiles1)
        sol.append(med_sol)
        interp_func_list.append(interp_func_dict)
        fig1.savefig(f'Files_n_figures/UVB_new_Q{q}_{file_name}_OVI.png')

        plt.close()

        # print('Including O VI \n')
        # fig2=corner.corner(flat_samples2, labels=labels, quantiles=[0.16, 0.5, 0.84],show_titles=True, title_kwargs={"fontsize": 12},verbose=True)
        # fig2.savefig(f'Files_n_figures/UVB_Q{q}_inc_OVI.png')
        print('\n')

    return sol,interp_func_list


obs_col_den=array([observations[i][0]  for i in observations.keys()])
col_den_error=array([observations[i][1]  for i in observations.keys()])


def interp_func(q):

    func_dict={}
    
    hdu=fits.open(f'Data/UV_Q_new/component_III_nH_UV_Q{q}_col_density_param.fits')
    data=Table(hdu[1].data)

    mask=data['H']!=0
    data=data[mask]

    log_nH=data['log_nH']

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


def chi_sq(sol,interp_func_dict,dof,logZ_ref=-1,inc_OVI=True):

    log_nH=sol[0]
    log_Z=sol[3]

    model_col_den=array([interp_func_dict[i](log_nH)+log_Z-logZ_ref for i in observations.keys()])

    if inc_OVI==True:

        chi_sq_val=sum((((model_col_den-obs_col_den)/(col_den_error))**2))
    
    else:
        chi_sq_val=sum((((model_col_den-obs_col_den)/(col_den_error))**2)[:-1])

    return chi_sq_val/dof

interp_func_Q=[interp_func(q) for q in range(14,14+n)]

ions_to_use1=['Si+', 'Si+2','C+', 'C+2']
ions_to_use2=['Si+', 'Si+2','C+', 'C+2','O+5']

sol_exc_OVI , interp_func_list = run_UVB(ions_to_use1,'excluding')
sol_inc_OVI = run_UVB(ions_to_use2,'including')[0]

# sol_exc_OVI=[[-2.96, 0.06, 0.07, -0.9, 0.09, 0.1], [-3.05, 0.06, 0.07, -0.96, 0.09, 0.11], [-3.13, 0.07, 0.07, -0.99, 0.1, 0.11], [-3.21, 0.07, 0.07, -1.04, 0.1, 0.12], [-3.29, 0.07, 0.07, -1.07, 0.1, 0.12], [-3.38, 0.07, 0.07, -1.12, 0.1, 0.12], [-3.45, 0.07, 0.07, -1.15, 0.1, 0.12]]
# sol_inc_OVI=[[-3.62, 0.01, 0.01, -1.21, 0.04, 0.04], [-3.74, 0.01, 0.02, -1.3, 0.04, 0.04], [-3.86, 0.02, 0.02, -1.38, 0.04, 0.04], [-3.97, 0.02, 0.02, -1.46, 0.04, 0.04], [-4.05, 0.02, 0.02, -1.51, 0.04, 0.04], [-4.15, 0.02, 0.02, -1.59, 0.04, 0.04], [-4.22, 0.02, 0.02, -1.64, 0.04, 0.04]]


chi_sq_exc1=[chi_sq(sol_exc_OVI[i],interp_func_Q[i],dof=2) for i in range(n)]
chi_sq_inc1=[chi_sq(sol_inc_OVI[i],interp_func_Q[i],dof=3) for i in range(n)]

chi_sq_exc2=[chi_sq(sol_exc_OVI[i],interp_func_Q[i],dof=2,inc_OVI=False) for i in range(n)]
chi_sq_inc2=[chi_sq(sol_inc_OVI[i],interp_func_Q[i],dof=3,inc_OVI=False) for i in range(n)]

x=linspace(1,n,n)
x_ticks=[f'Q{int(i+13)}' for i in x]

plt.figure()

plt.subplot(2,2,1)

plt.title('Including OVI data point for '+r'$\chi^2$ calculation')
plt.scatter(x,chi_sq_inc1,label='including OVI',c='green')
plt.plot(x,chi_sq_inc1,c='green',ls='--')
plt.legend()
plt.ylabel(r'$\chi^2$',labelpad=15)
plt.xticks(x,x_ticks)

plt.subplot(2,2,2)

plt.title('Including OVI data point for '+r'$\chi^2$ calculation')
plt.scatter(x,chi_sq_exc1,label='excluding OVI',c='green')
plt.plot(x,chi_sq_exc1,ls='--',c='green')
plt.legend()
plt.xticks(x,x_ticks)

plt.subplot(2,2,3)

plt.title('Excluding OVI data point for '+r'$\chi^2$ calculation')
plt.scatter(x,chi_sq_inc2,label='including OVI',c='green')
plt.plot(x,chi_sq_inc2,c='green',ls='--')
plt.legend()
plt.xlabel('Q',labelpad=15)
plt.ylabel(r'$\chi^2$',labelpad=15)
plt.xticks(x,x_ticks)

plt.subplot(2,2,4)

plt.title('Excluding OVI data point for '+r'$\chi^2$ calculation')
plt.scatter(x,chi_sq_exc2,label='excluding OVI',c='green')
plt.plot(x,chi_sq_exc2,ls='--',c='green')
plt.legend()
plt.xlabel('Q',labelpad=15)
plt.xticks(x,x_ticks)



nH_exc=[]
nH_l_err_exc=[]
nH_h_err_exc=[]
Z_exc=[]
Z_l_err_exc=[]
Z_h_err_exc=[]

nH_inc=[]
nH_l_err_inc=[]
nH_h_err_inc=[]
Z_inc=[]
Z_l_err_inc=[]
Z_h_err_inc=[]

for i in range(n):

    nH_exc.append(sol_exc_OVI[i][0])
    nH_l_err_exc.append(sol_exc_OVI[i][1])
    nH_h_err_exc.append(sol_exc_OVI[i][2])
    Z_exc.append(sol_exc_OVI[i][3])
    Z_l_err_exc.append(sol_exc_OVI[i][4])
    Z_h_err_exc.append(sol_exc_OVI[i][5])

    nH_inc.append(sol_inc_OVI[i][0])
    nH_l_err_inc.append(sol_inc_OVI[i][1])
    nH_h_err_inc.append(sol_inc_OVI[i][2])
    Z_inc.append(sol_inc_OVI[i][3])
    Z_l_err_inc.append(sol_inc_OVI[i][4])
    Z_h_err_inc.append(sol_inc_OVI[i][5])


x=linspace(1,n,n)
x_ticks=[f'Q{int(i+13)}' for i in x]

plt.figure()
plt.title('logN values for different UV backgrounds')
plt.errorbar(x,nH_exc,yerr=(nH_l_err_exc,nH_h_err_exc),capsize=3,label='Excluding O VI',fmt='o',c='red')
plt.errorbar(x,nH_inc,yerr=(nH_l_err_inc,nH_h_err_inc),capsize=3,label='Including O VI',fmt='o',c='green')
plt.plot(x,nH_exc,ls='--',c='red')
plt.plot(x,nH_inc,ls='--',c='green')
plt.xticks(x,x_ticks)
plt.xlabel('Q value')
plt.ylabel('log H')
plt.legend()

plt.figure()
plt.title('logZ values for different UV backgrounds')

plt.errorbar(x,Z_exc,yerr=(Z_l_err_exc,Z_h_err_exc),capsize=3,label='Excluding O VI',fmt='o',c='red')
plt.errorbar(x,Z_inc,yerr=(Z_l_err_inc,Z_h_err_inc),capsize=3,label='Including O VI',fmt='o',c='green')
plt.plot(x,Z_exc,ls='--',c='red')
plt.plot(x,Z_inc,ls='--',c='green')
plt.xticks(x,x_ticks)
plt.xlabel('Q value')
plt.ylabel('log Z')
plt.legend()

markers=['o','s','+','x','D','^','p']

plt.figure()

for i in range(n):

    plt.errorbar(nH_exc[i],Z_exc[i],yerr=([Z_l_err_exc[i]],[Z_h_err_exc[i]]),xerr=([nH_l_err_exc[i]],[nH_h_err_exc[i]]),c='red',marker=markers[i],label=f'Q{14+i}',capsize=3)
    plt.scatter(nH_exc[i],Z_exc[i],c='red',marker=markers[i],s=100)

plt.plot(nH_exc,Z_exc,c='red',ls='--')

plt.title('Cloudy model solutions (excluding OVI)')
plt.xlabel('log H')
plt.ylabel('log Z')
plt.legend()

plt.figure()

for i in range(n):

    plt.errorbar(nH_inc[i],Z_inc[i],yerr=([Z_l_err_inc[i]],[Z_h_err_inc[i]]),xerr=([nH_l_err_inc[i]],[nH_h_err_inc[i]]),c='green',marker=markers[i],label=f'Q{14+i}',capsize=3)
    plt.scatter(nH_inc[i],Z_inc[i],c='green',marker=markers[i],s=100)

plt.plot(nH_inc,Z_inc,c='green',ls='--')

plt.title('Cloudy model solutions (including OVI)')
plt.xlabel('log H')
plt.ylabel('log Z')
plt.legend()

plt.show()

quit()

def plot_UVB():

    sol={14:sol_14, 15:sol_15, 16:sol_16, 17:sol_17, 18:sol_18, 19:sol_19, 20:sol_20}

    ions=['Si+2', 'Si+3','C+2', 'C+3','O+6']

    # observations={'Si+2':[12.53,0.07], 'Si+3':[12.69,0.06],'C+2':[13.49,0.05], 'C+3':[13.65,0.03],'O+6':[13.84,0.03]}
    observations={'Si+2':[12.53,0.21], 'Si+3':[12.69,0.19],'C+2':[13.49,0.21], 'C+3':[13.65,0.23],'O+6':[13.84,0.22]}

    obs_col_den=[observations[i][0]  for i in ions]
    col_den_error=[observations[i][1]  for i in ions]

    col_den_exc_OVI=[]
    col_den_inc_OVI=[]

    for q in range(14,21):

        # print(f'UV_q = {q} \n')

        hdu=fits.open(f'Data/component_III_1d_Q{q}_col_density_param.fits')
        data=Table(hdu[1].data)

        mask=data['H']!=0
        data=data[mask]

        log_nH=data['log_nH']

        logZ_ref=-1

        lognH_range=[-5,2]
        logZ_range=[-3,2]

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
                f=interp1d(log_nH,log_col_den,fill_value='extrapolate')
                func_dict[i]=f

            return func_dict

        interp_func_dict=interp_func()

        nH1,Z1=sol[q].sol_exc_OVI
        nH2,Z2=sol[q].sol_inc_OVI

        print(nH1,Z1,'\n',nH2,Z2, '\n')

        median_col_den_exc_OVI=[interp_func_dict[i](nH1)+Z1-logZ_ref for i in observations.keys()]
        median_col_den_inc_OVI=[interp_func_dict[i](nH2)+Z2-logZ_ref for i in observations.keys()]

        col_den_exc_OVI.append(median_col_den_exc_OVI)
        col_den_inc_OVI.append(median_col_den_inc_OVI)

    ions_all=[f'Si {toRoman(2)}', f'Si {toRoman(3)}',f'C {toRoman(2)}', f'C {toRoman(3)}',f'O {toRoman(6)}']
    x=linspace(1,len(ions_all),len(ions_all))

    plt.figure()

    plt.errorbar(x,obs_col_den,c='red',yerr=col_den_error, fmt='o',capsize=3,label='Observed')

    for i in range(14,21):
        plt.plot(x,col_den_exc_OVI[i-14],ls='--',label=f'Q={i}')
    
    plt.xticks(x,ions_all,fontsize=20)
    plt.yticks(arange(12.5,14.1,0.5),fontsize=20)
    plt.ylabel(r'$\mathbf{log \ (N \ {cm}^{-2})}$',labelpad=15)
    plt.xlabel(r'$\mathbf{Ions}$',labelpad=15)
    plt.legend()
    plt.title('Solutions for different UV backgrounds (Excluding O VI)')
    
    plt.figure()

    plt.errorbar(x,obs_col_den,c='red',yerr=col_den_error, fmt='o',capsize=3,label='Observed')

    for i in range(14,21):
        plt.plot(x,col_den_inc_OVI[i-14],ls='-.',label=f'Q={i}')

    plt.xticks(x,ions_all,fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylabel(r'$\mathbf{log \ (N \ {cm}^{-2})}$',labelpad=15)
    plt.xlabel(r'$\mathbf{Ions}$',labelpad=15)
    plt.legend()
    plt.title('Solutions for different UV backgrounds (Including O VI)')


plot_UVB()
plt.show()




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
        int_col_den1.append(f(sample1[0])+sample1[1]-logZ_ref)
        int_col_den2.append(f(sample2[0])+sample2[1]-logZ_ref)

    # plt.plot(x,int_col_den1,c='orange',alpha=0.1)
    if i!=inds[-1]:
        plt.plot(x,int_col_den1,c='orange',alpha=0.1)
        plt.plot(x,int_col_den2,c='green',alpha=0.1)
    
    else:
        plt.plot(x,int_col_den1,c='orange',alpha=0.4,label='Model samples (excluding O VI)')
        plt.plot(x,int_col_den2,c='green',alpha=0.4,label='Model samples (including O VI)')


median_col_den_exc_OVI=[interp_func_dict[i](-3.54)-1.37-logZ_ref for i in observations.keys()]
median_col_den_inc_OVI=[interp_func_dict[i](-4.07)-1.25-logZ_ref for i in observations.keys()]

obs_col_den=[observations[i][0]  for i in ions]
col_den_error=[observations[i][1]  for i in ions]


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


























# class UV_sol():

#     def __init__(self,nH1,Z1,nH2,Z2):

#         self.sol_exc_OVI=[nH1[1][1],Z1[1][1]]
#         self.err_exc_OVI=[[nH1[1][1]-nH1[0][1],nH1[2][1]-nH1[1][1]],[Z1[1][1]-Z1[0][1],Z1[2][1]-Z1[1][1]]]
        
#         self.sol_inc_OVI=[nH2[1][1],Z2[1][1]]
#         self.err_inc_OVI=[[nH2[1][1]-nH2[0][1],nH2[2][1]-nH2[1][1]],[Z2[1][1]-Z2[0][1],Z2[2][1]-Z2[1][1]]]


# sol_14=UV_sol([(0.16, -3.1650557555497354), (0.5, -3.007405534604607), (0.84, -2.816370144285358)],
# [(0.16, -1.2140547198754161), (0.5, -1.0495140997870898), (0.84, -0.865016578770266)],
# [(0.16, -3.5458092604638782), (0.5, -3.481049618177483), (0.84, -3.4174008591564937)],
# [(0.16, -0.981496687934417), (0.5, -0.8482992055550067), (0.84, -0.7141726932442224)])

# sol_15=UV_sol([(0.16, -3.264071557288731), (0.5, -3.0953970136147806), (0.84, -2.907400979167775)],
# [(0.16, -1.2711167941633967), (0.5, -1.1049544948109422), (0.84, -0.922653056911785)],
# [(0.16, -3.6653150254647415), (0.5, -3.6013363754628656), (0.84, -3.532952290308554)],
# [(0.16, -1.0505940862608458), (0.5, -0.9196313738607506), (0.84, -0.7858551059587006)])

# sol_16=UV_sol([(0.16, -3.3586884710721234), (0.5, -3.1897035732214984), (0.84, -2.9949410809283017)],
# [(0.16, -1.3361192253598315), (0.5, -1.1681537935350923), (0.84, -0.9834892023136539)],
# [(0.16, -3.779785622116329), (0.5, -3.7121759902296674), (0.84, -3.6429655453795178)],
# [(0.16, -1.125879718136662), (0.5, -0.99577248219591), (0.84, -0.8599572623989533)])

# sol_17=UV_sol([(0.16, -3.460224639748311), (0.5, -3.2859923202343104), (0.84, -3.0937758285667165)],
# [(0.16, -1.3927537317687266), (0.5, -1.2234736226282785), (0.84, -1.039727367275296)],
# [(0.16, -3.8847712132961356), (0.5, -3.816254080753127), (0.84, -3.7464886192792592)],
# [(0.16, -1.193564583444682), (0.5, -1.0590634870824858), (0.84, -0.9278276962239561)])

# sol_18=UV_sol([(0.16, -3.5413240223361955), (0.5, -3.3664139941772593), (0.84, -3.171298685408941)],
# [(0.16, -1.4376091762330279), (0.5, -1.270556849562071), (0.84, -1.086166248548544)],
# [(0.16, -3.9721856088616514), (0.5, -3.9007531788614678), (0.84, -3.828194675504072)],
# [(0.16, -1.248828168014952), (0.5, -1.1179046710509706), (0.84, -0.9831465839294989)])

# sol_19=UV_sol([(0.16, -3.6327152974875823), (0.5, -3.454752197385604), (0.84, -3.2582490368087234)],
# [(0.16, -1.5021454644310077), (0.5, -1.3322621439188398), (0.84, -1.1417100780903657)],
# [(0.16, -4.066173359313065), (0.5, -3.9942821834130458), (0.84, -3.9205176438637923)],
# [(0.16, -1.3257512243380025), (0.5, -1.1924538569991192), (0.84, -1.058315432865459)])

# sol_20=UV_sol([(0.16, -3.713196981676712), (0.5, -3.5313040472917514), (0.84, -3.3300844958145666)],
# [(0.16, -1.5410646153028928), (0.5, -1.3721606948585217), (0.84, -1.1832223449169355)],
# [(0.16, -4.139886829340733), (0.5, -4.065465366851797), (0.84, -3.9895003157514157)],
# [(0.16, -1.3842845621923168), (0.5, -1.2517900855782296), (0.84, -1.1181521168276485)])









