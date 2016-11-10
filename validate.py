import aimflex
import corner
import numpy as np
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator
from astropy.io import fits
from astropy.table import Table,hstack
from os import path
import datetime
import sys

if len(sys.argv)<2:
	outdir = '../../flexion_development/val_out'
else:
	outdir = sys.argv[1]


field_data_raw = Table.read('true_fields.txt',format='ascii.fixed_width')
e_sigma = 0.2		# normal
alpha_rng = (2,15)	# uniform
logI_rng = (1,2)	# uniform
index_rng=(0.25,5)	# uniform in log
noise_sig = 1.		# normal

sample_size=2
break_size =100

dims = (101,101)
x,y =np.meshgrid(np.linspace(-(dims[0]-1)/2,(dims[0]-1)/2,dims[0]),
                 np.linspace(-(dims[1]-1)/2,(dims[1]-1)/2,dims[1]))
p=None
model = aimflex.get_AIM()

E = np.random.randn(sample_size,2)*e_sigma
A = np.random.uniform(*alpha_rng,size=sample_size)
logI = np.random.uniform(*logI_rng,size=sample_size)
index = np.power(10,np.random.uniform(np.log10(index_rng[0]),
									  np.log10(index_rng[1]),
									  sample_size))

sample = np.sort(np.random.randint(0,len(field_data_raw),sample_size))


input_pars = Table(names=model.param_names)
best_pars = Table(names=[pn+'_best'for pn in model.param_names])
mean_pars = Table(names=[pn+'_mean'for pn in model.param_names])
median_pars = Table(names=[pn+'_median'for pn in model.param_names])
obj_tbl = Table(names=['objnum','sample_point'],dtype=('i4','i4',))

upper_error = Table(names=[pn+'_err_up'for pn in model.param_names])
lower_error = Table(names=[pn+'_err_lo'for pn in model.param_names])
sigma_error = Table(names=[pn+'_stdev'for pn in model.param_names])

lsq = Table(names=['lsq_value'])

for i,s in enumerate(sample):
	print i,datetime.datetime.now()
	if (i % break_size) == 0 and i>1:
		input_data = hstack([obj_tbl,input_pars])
		output_data = hstack([obj_tbl,
							  best_pars,mean_pars,median_pars,
							  sigma_error,upper_error,lower_error,
							  lsq])

		input_data.write(path.join(outdir,"input_data_{:03d}.txt".format(i-break_size)),format='ascii.fixed_width')
		output_data.write(path.join(outdir,"output_data_{}.txt"),format='ascii.fixed_width')


	model.c1_0 = 0.
	model.c2_0 = 0.

	model.E1_2 = E[i,0]
	model.E2_2 = E[i,1]
	model.alpha_2 = A[i]
	model.logI_2 = logI[i]
	model.index_2 = index[i]
	
	model.g1_0 = field_data_raw['g1'][s]
	model.g2_0 = field_data_raw['g2'][s]
	model.F1_0 = field_data_raw['F1'][s]
	model.F2_0 = field_data_raw['F2'][s]
	model.G1_0 = field_data_raw['G1'][s]
	model.G2_0 = field_data_raw['G2'][s]
	
	true_params = np.copy(model.parameters)
	
	true = model(x,y,p)
	data = true + aimflex.window_image(np.random.randn(*(true.shape)))
	weights = np.ones(true.shape)

	aimflex.set_gaussian_pars(data,model,weights=weights)
	aimflex.set_limits(data,model)

	### FITTING
	###
	samples, chisqs = aimflex.fit_image(model,data,weights,p,verbose=False)
	samples2D=samples.reshape((-1, model.parameters.size))

	# corner figure
	cfig = corner.corner(samples2D, labels=model.param_names,
					truths=true_params)
	cfig.savefig(path.join(outdir,"obj_{:03d}_corner.png".format(i)))
	cfig.clf()
	
	# results to save
	stats = aimflex.avg_err_corr(samples2D)
	best = np.argmin(chisqs)
	
	fit_pars = samples[best]
	chisq = chisqs[best]
	
	model.parameters = fit_pars
	
	outimg = model(x,y,p)

	# Image plots
	f, axarr = pl.subplots(2,2)
	axarr[0,0].imshow(data,cmap='bone')
	axarr[0,1].imshow(true,cmap='bone')
	axarr[1,0].imshow(outimg,cmap='bone')
	axarr[1,1].imshow(data-outimg,cmap='bone')

	axarr[0,0].set_title('Data')
	axarr[0,1].set_title('Input')
	axarr[1,0].set_title('Output')
	axarr[1,1].set_title('Residual')

	f.savefig(path.join(outdir,"obj_{:03d}_images.png".format(i)))
	f.clf()
	
	# Parameter track plots
	fig, axes = pl.subplots(model.parameters.size, 1, sharex=True, figsize=(8, 20))

	for k in range(model.parameters.size):
		axes[k].plot(samples[:, k].T, color="k", alpha=0.4)
		axes[k].yaxis.set_major_locator(MaxNLocator(5))
		axes[k].axhline(true_params[k], color="b", lw=2)
		axes[k].set_ylabel(model.param_names[k])

	fig.tight_layout(h_pad=0.0)
	fig.savefig(path.join(outdir,"obj_{:03d}_pars.png".format(i)))
	fig.clf()
	
	del fig, cfig, f

	# Output
	obj_tbl.add_row([i,s])
	
	input_pars.add_row(true_params)
	best_pars.add_row(fit_pars)
	mean_pars.add_row(stats[0])
	sigma_error.add_row(stats[1])

	median_pars.add_row(stats[2][1])
	upper_error.add_row(stats[2][2] - stats[2][1])
	lower_error.add_row(stats[2][2] - stats[2][1])

	lsq.add_row(chisq)


# 	input_data.add_row(true_params)
# 	output_data.add_row(mcmc_results[:,0])
# 	chisq.add_row([aimflex.leastsquare(data, model, weights, x, y, p)])
# 	samp_tbl.add_row([s])
# 	obj_tbl.add_row([i])
# 	upper_error.add_row(mcmc_results[:,1])
# 	lower_error.add_row(mcmc_results[:,1])
	
	pl.close("all")


input_data.write(path.join(outdir,"input_data_master.txt"),format='ascii.fixed_width')
output_data.write(path.join(outdir,"output_data_master.txt"),format='ascii.fixed_width')
print 'Finished!!',datetime.datetime.now()

