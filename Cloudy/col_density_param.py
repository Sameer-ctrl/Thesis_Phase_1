from numpy import loadtxt
from astropy.io import ascii

grid=loadtxt('../Cloudy_runs/inference1d_grid.txt',dtype=str)
col_density=ascii.read('../Cloudy_runs/inference1d_column_density.txt')


nH=grid[:,6].astype(float)
# Z=grid[:,7].astype(float)
# T=grid[:,8].astype(float)

# params=tuple(zip(nH,Z,T))
params=nH

col_density.add_column(params,name='parameters')
col_density.add_column(nH,name='nH')
# col_density.add_column(Z,name='Z')
# col_density.add_column(T,name='T')

col_density.write('../Cloudy_runs/Inference1d_col_density_param.fits',overwrite=True)