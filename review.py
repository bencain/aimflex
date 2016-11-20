import aimflex
import numpy as np
import matplotlib.pyplot as pl
from astropy.table import Table,hstack
from os import path

datadir='../../flexion_development/val_out_20161116/'

true=Table.read(path.join(datadir,'input_data_master.txt'),format='ascii.fixed_width')
fit=Table.read(path.join(datadir,'output_data_master.txt'),format='ascii.fixed_width')


# nc = np.ceil(np.sqrt(len(true.colnames[2:]))).astype(int)
# fig, axes = pl.subplots(nc,nc)


# Focus on the shape parameters
shapep = ['E1_2','g1_0','E2_2','g2_0',
		  'F1_0','F2_0','G1_0','G2_0']
fig, axes = pl.subplots(2,4, figsize=(12, 6))
nc=4

ns=100
samp = np.random.randint(0,true[shapep[0]].size,ns)
n=np.arange(ns)
samp=np.arange(true[shapep[0]].size)


lim=[1,1,1,1,0.1,0.1,0.1,0.1]
# for k,p in enumerate(true.colnames[2:]):
for k,p in enumerate(shapep):
	tmin = min(true[p][samp])
	tmax = max(true[p][samp])
	axes[k//nc,k%nc].errorbar(true[p][samp],fit[p+'_best'][samp],yerr=[fit[p+'_err_lo'][samp],fit[p+'_err_up'][samp]],linestyle='none',fmt='b.')
	axes[k//nc,k%nc].set_xlabel(p+' input')
	axes[k//nc,k%nc].set_ylabel(p+' output')
	axes[k//nc,k%nc].plot([tmin,tmax],[tmin,tmax],'k')
	axes[k//nc,k%nc].set_xlim((lim[k],-lim[k]))
	axes[k//nc,k%nc].set_ylim((lim[k],-lim[k]))

fig.savefig(path.join(datadir,"shape_input_vs_output.png"))

pl.clf()

# Plot true vs 
#	best
#	mean
#	median
#	mode
# for each shape parameter



for i,ip in enumerate(shapep):
	fig, axes = pl.subplots(2,2, figsize=(8, 8))
	for k,tag in enumerate(['_best','_mean','_median','_mode']):
		op = ip+tag
# 		axes[k//2,k%2].errorbar(n,fit[op][samp]-true[ip][samp],yerr=[fit[ip+'_err_lo'][samp],fit[ip+'_err_up'][samp]],linestyle='none',fmt='b.')	
# 		axes[k//2,k%2].errorbar(n,fit[op][samp]-true[ip][samp],yerr=fit[ip+'_stdev'][samp],linestyle='none',fmt='r.')	

# 		axes[k//2,k%2].scatter(fit[ip+'_stdev'][samp],fit[op][samp]-true[ip][samp],marker='.')	
# 		axes[k//2,k%2].set_xlabel(ip+'_stdev')

		axes[k//2,k%2].scatter(fit['alpha_2'+tag][samp],(fit[op][samp]-true[ip][samp])/fit[ip+'_stdev'][samp],marker='x')
		axes[k//2,k%2].set_ylabel(op+' - '+ip+'/sigma')

		axes[k//2,k%2].set_ylim((lim[i],-lim[i]))

# 	pl.show()
# 	pl.clf()
	fig.savefig(path.join(datadir,"shape_metric_{}.png".format(ip)))
	pl.clf()

		
