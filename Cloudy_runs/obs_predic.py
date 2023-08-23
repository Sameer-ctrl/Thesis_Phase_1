import matplotlib.pyplot as plt
from numpy import *
from roman import toRoman


plt.style.use('../Python/my_style.mpl')


ions=[f'Si {toRoman(2)}', f'Si {toRoman(3)}',f'C {toRoman(2)}', f'C {toRoman(3)}', f'O {toRoman(6)}']
x=linspace(1,len(ions),len(ions))


obs_col_density=array([12.53,12.69,13.49,13.65,13.84])
col_density_error=array([0.07,0.06,0.05,0.03,0.03])

predicted_all_ions=log10([4.43676e+11, 6.40882e+11, 1.18221e+14, 1.23953e+14, 4.68142e+13])
predicted_without_OVI=log10([3.29517e+12, 2.07164e+12,1.02575e+14, 3.75366e+13, 4.88502e+11])

fig=plt.figure(figsize=(11,11))

plt.errorbar(x,obs_col_density,c='red',yerr=col_density_error, fmt='o',capsize=3,label='Observed')
plt.plot(x,predicted_all_ions,label=r'All ions $(T={10}^{4.49}\ K)$',ls='--')
plt.plot(x,predicted_without_OVI,label=r'Excluding OVI $(T={10}^{4.38}\ K)$',ls='--')
plt.xticks(x,ions,fontsize=20)
plt.yticks(fontsize=20)
plt.ylabel(r'$\mathbf{log \ (N \ {cm}^{-2})}$',labelpad=15)
plt.xlabel(r'$\mathbf{Ions}$',labelpad=15)
plt.legend(loc='upper left')
plt.savefig('Files_n_figures/Observed_and_predicted.png')
plt.show()
