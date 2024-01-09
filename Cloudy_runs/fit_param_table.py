from numpy import *
import os

qso='pg1121'
z_abs=0.192393

file=f'../VPfit/{qso}/z={z_abs}/fit_params.txt'

qso_list=loadtxt('../Python/Data/qso_list.txt',dtype=str)
qso_dict=dict(zip(qso_list[:,1],qso_list[:,0]))
qso_label=qso_dict[qso]


with open(file) as f:
    data=f.read()
    data=data.replace('A','')
    print(data)
    print('\n')

with open('temp_param_file.txt','w') as f:
    f.write(data)

param_file=loadtxt('temp_param_file.txt',dtype=str,comments=('>','#'))

ions=param_file[:,0]
z=param_file[:,1].astype(float)
z_err=param_file[:,2].astype(float)
b=param_file[:,3].astype(float)
b_err=param_file[:,4].astype(float)
logN=param_file[:,5].astype(float)
logN_err=param_file[:,6].astype(float)


v_z=lambda z : 3e5*(((1+round(z,6))**2-1)/((1+round(z,6))**2+1))
err_vz=lambda z,z_err: 4*3e5*((1+round(z,6))/(((1+round(z,6))**2)+1)**2)*round(z_err,6)

def ion_format(ion):

    for i,s in enumerate(ion):
        if s.isupper():
            if i==0:
                pass
            else:
                break

    atom=ion[:i]
    n=ion[i:]

    return f'\ion{{{atom}}}{{{n.lower()}}}'

v_abs=v_z(z_abs)
v=zeros(len(ions))
v_err=zeros(len(ions))

for i in range(len(ions)):

    v[i]=round(v_z(z[i])-v_abs)
    v_err[i]=round(err_vz(z[i],z_err[i]))


print('\\newpage\n\n\\begin{landscape}\n\n\\begin{figure}\n    \centering\n    \\vspace{-20mm}\n    \hspace*{-35mm}\n    \captionsetup{oneside,margin={0cm,35mm}}')
print(f'    \includegraphics[width=1.25\linewidth]{{System-Plots/{qso_label}_z={z_abs}_sys_plot.png}}\n\end{{figure}}\n\n\end{{landscape}}\n\n')
print('\\begin{center} \n\n\\begin{tabular}{cccc} \n\n    \hline \hline \\tabularnewline \n    \head{Ion} & \head{v (km s\\textsuperscript{$\mathbf{-1}$})} & \head{b (km s\\textsuperscript{$\mathbf{-1}$})} & \head{log [N cm\\textsuperscript{$\mathbf{-2}$}]}\n    \\tabularnewline \\tabularnewline \hline \\tabularnewline \n ')

        
for i in range(len(ions)):
    print(f'    {ion_format(ions[i])}   &    {v[i]} $\pm$ {v_err[i]}   &    {round(b[i])} $\pm$ {round(b_err[i])}    &     {round(logN[i],2)} $\pm$ {round(logN_err[i],2)} \\\\')

print('\n    \\tabularnewline \hline \hline \\tabularnewline \n\n\end{tabular}\n\n\end{center}')



    










os.remove('temp_param_file.txt')