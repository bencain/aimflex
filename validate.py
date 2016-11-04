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

sample_size=100

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

sample = np.sort(np.random.random_integers(0,len(field_data_raw)-1,
											sample_size))


input_data = Table(names=model.param_names)
output_data = Table(names=model.param_names)
chisq = Table(names=['chisq'])
samp_tbl = Table(names=['sample_point'],dtype=('i4',))
obj_tbl = Table(names=['objnum'],dtype=('i4',))
upper_error = Table(names=[pn+'_err_up'for pn in model.param_names])
lower_error = Table(names=[pn+'_err_lo'for pn in model.param_names])

# pl.ion()
for i,s in enumerate(sample):
	print i,datetime.datetime.now()

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
	samples = aimflex.fit_image(model,data,weights,p,verbose=False)
# 	samples = np.random.randn(500,13)

	cfig = corner.corner(samples, labels=model.param_names,
					truths=true_params)
	cfig.savefig(path.join(outdir,"obj_{:03d}_corner.png".format(i)))
	cfig.clf()
	
	
	# tuples of (median, dx+, dx-) for 68% confidence contours
	conf = 0.68
	mcmc_results = np.array(
					map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
					zip(*np.percentile(samples, [50-conf*100/2, 
												 50, 
												 50+conf*100/2],axis=0)))
					)

	model.parameters = mcmc_results[:,0]
	outimg = model(x,y,p)

	f, axarr = pl.subplots(2,2)
	axarr[0,0].imshow(data)
	axarr[0,1].imshow(true)
	axarr[1,0].imshow(outimg)
	axarr[1,1].imshow(data-outimg)

	axarr[0,0].set_title('Data')
	axarr[0,1].set_title('Input')
	axarr[1,0].set_title('Output')
	axarr[1,1].set_title('Residual')

	f.savefig(path.join(outdir,"obj_{:03d}_images.png".format(i)))
	f.clf()
	
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

	input_data.add_row(true_params)
	output_data.add_row(mcmc_results[:,0])
	chisq.add_row([aimflex.leastsquare(data, model, weights, x, y, p)])
	samp_tbl.add_row([s])
	obj_tbl.add_row([i])
	upper_error.add_row(mcmc_results[:,1])
	lower_error.add_row(mcmc_results[:,1])
	
	pl.close("all")

print 'Finished!!',datetime.datetime.now()

# 	pl.imshow(data)
# 	print model
# 	raw_input("Next?")

input_data = hstack([obj_tbl,input_data,samp_tbl])
output_data = hstack([obj_tbl,output_data,chisq,upper_error,lower_error,samp_tbl])

input_data.write(path.join(outdir,"input_data.txt"),format='ascii.fixed_width')
output_data.write(path.join(outdir,"output_data.txt"),format='ascii.fixed_width')
