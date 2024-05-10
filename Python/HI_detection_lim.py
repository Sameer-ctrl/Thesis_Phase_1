from numpy import *
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.constants import * 
from scipy.integrate import quad


rest_wave=1215.671
wave=1215.671e-10   # Ang

f=0.4164
gamma=3.363e34
tau=0.1

N=logspace(11,14,100)

b=(sqrt(pi)*(elementary_charge**2)/(electron_mass*speed_of_light))*wave*f*(10**N)*(1/tau)

plt.plot(N,b)
plt.show()
quit()

# def N(b):

#     a=(wave*gamma)/(4*pi*b)
#     u=(-speed_of_light/b)*(1-(wave/rest_wave))

#     def func(x):

#         return (a/pi)*(exp(-x**2)/((u-x)**2+(a**2)))

#     H=quad(func,-inf,inf)[0]

#     col_den=(tau*b)/(0.01498*f*wave*H)

#     return col_den

# b=linspace(10000,200000,500)

# N=[N(x) for x in b]

# plt.plot(log10(N),log10(b))

# plt.show()