import aimflex
import numpy as np
import matplotlib.pyplot as pl
from astropy.table import Table,hstack
from os import path

datadir='../../flexion_development/val_out_20161103/'

true=Table.read(path.join(datadir,'input_data.txt'),format='ascii.fixed_width')
fit=Table.read(path.join(datadir,'output_data.txt'),format='ascii.fixed_width')

nc = len(true.colnames)//2 + len(true.colnames)%2

fig, axes = pl.subplots(2,nc)

for k,p in enumerate(true.colnames):
	tmin = min(true[p])
	tmax = max(true[p])
	axes[k//nc,k%nc].errorbar(true[p],fit[p],yerr=[fit[p+'_err_lo'],fit[p+'_err_up']],linestyle='none',fmt='b.')
	axes[k//nc,k%nc].set_xlabel(p+' input')
	axes[k//nc,k%nc].set_ylabel(p+' output')
	axes[k//nc,k%nc].plot([tmin,tmax],[tmin,tmax],'k')

# fig.tight_layout(h_pad=0.0)
# fig.savefig(path.join(outdir,"obj_{:03d}_pars.png".format(i)))
# fig.clf()
pl.show()
