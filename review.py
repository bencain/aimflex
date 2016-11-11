import aimflex
import numpy as np
import matplotlib.pyplot as pl
from astropy.table import Table,hstack
from os import path

datadir='../../flexion_development/val_out_20161104/'
datadir='/Users/bcain/Desktop/tmp/'

true=Table.read(path.join(datadir,'input_data_210.txt'),format='ascii.fixed_width')
fit=Table.read(path.join(datadir,'output_data_210.txt'),format='ascii.fixed_width')


# nc = np.ceil(np.sqrt(len(true.colnames[2:]))).astype(int)
# fig, axes = pl.subplots(nc,nc)



shapep = ['E1_2','g1_0','E2_2','g2_0',
		  'F1_0','F2_0','G1_0','G2_0']
fig, axes = pl.subplots(2,4)
nc=4

# for k,p in enumerate(true.colnames[2:]):
for k,p in enumerate(shapep):
	tmin = min(true[p])
	tmax = max(true[p])
	axes[k//nc,k%nc].errorbar(true[p],fit[p+'_best'],yerr=[fit[p+'_err_lo'],fit[p+'_err_up']],linestyle='none',fmt='b.')
	axes[k//nc,k%nc].set_xlabel(p+' input')
	axes[k//nc,k%nc].set_ylabel(p+' output')
	axes[k//nc,k%nc].plot([tmin,tmax],[tmin,tmax],'k')

# fig.tight_layout(h_pad=0.0)
# fig.savefig(path.join(outdir,"obj_{:03d}_pars.png".format(i)))
# fig.clf()
pl.show()
