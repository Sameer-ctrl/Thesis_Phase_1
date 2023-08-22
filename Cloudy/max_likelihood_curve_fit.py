from numpy import *
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit



x=linspace(0,1,20)
y=sin(1.5*x**2)+(4.57*(x**3))


noise=random.normal(0,0.25,len(x))
sig=mean(abs(noise))

y_noise=y+noise


def fit_func(x,a,b):

    return sin(a*x**2)+(b*(x**3))


def chi_sq(param):

    return sum((y_noise-fit_func(x,param[0],param[1]))**2)/(2*(sig**2))

def like(chi_sq):

    return e**(-chi_sq)

def next_guess(current_param):

    r=0.01
    new_param1=random.uniform(current_param[0]-r,current_param[0]+r)
    new_param2=random.uniform(current_param[1]-r,current_param[1]+r)

    return [new_param1,new_param2]

def like_ratio(current_param,proposed_param):

    return e**(chi_sq(current_param)-chi_sq(proposed_param))

a0=1
b0=1
n=10000

chi_sq_init=chi_sq([a0,b0])
current_param=[a0,b0]
a=[]
b=[]
likelihood=[]


for i in range(n):
    a.append(current_param[0])
    b.append(current_param[1])
    likelihood.append(like(chi_sq(current_param)))
    new_param=next_guess(current_param)
    l_ratio=like_ratio(current_param,new_param)

    r=random.random(1)[0]

    if l_ratio>=r:
        current_param=new_param
    
    else:
        pass

fit=curve_fit(fit_func,x,y_noise)
i=argmax(likelihood)
print(f'MCMC: a = {a[i]} b = {b[i]}')
print(f'chi_sq: a = {fit[0][0]} b = {fit[0][1]}')


plt.subplot(2,2,1)
plt.hist(a)

plt.subplot(2,2,3)
plt.scatter(a,b)

plt.subplot(2,2,4)
plt.hist(b)

plt.show()







# plt.plot(x,y,c='green')
# plt.plot(x,fit_func(x,a[i],b[i]),ls='--',c='red',label='MCMC')
# plt.plot(x,fit_func(x,fit[0][0],fit[0][1]),ls='--',c='orange',label='chi_sq')
# plt.scatter(x,y_noise)
# plt.legend()
plt.show()




