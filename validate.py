import aimflex
import corner
import sep
import datetime
import sys

import numpy as np
import matplotlib.pyplot as pl

from matplotlib.ticker import MaxNLocator
from astropy.io import fits
from astropy.table import Table,hstack
from os import path

if len(sys.argv)<2:
	outdir = '../../flexion_development/val_out'
else:
	outdir = sys.argv[1]


field_data_raw = Table.read('true_fields.txt',format='ascii.fixed_width')
e_sigma = 0.2		# normal
alpha_rng = (5,15)	# uniform
logI_rng = (0.5,2)	# uniform
index_rng=(0.2,5)	# uniform
noise_sig = 1.		# normal

cut_factor=3.0		# multiples of measured, lensed semimajor axis size

sample_size=100
break_size =10

dims = (201,201)
x,y =np.meshgrid(np.linspace(-(dims[0]-1)/2,(dims[0]-1)/2,dims[0]),
                 np.linspace(-(dims[1]-1)/2,(dims[1]-1)/2,dims[1]))
p=None
inmodel = aimflex.get_aim(flipped=True)
outmodel= aimflex.get_aim(gaussian=True)
parmask = [pn in outmodel.param_names for pn in inmodel.param_names]

E = np.random.randn(sample_size,2)*e_sigma
A = np.random.uniform(*alpha_rng,size=sample_size)
logI = np.random.uniform(*logI_rng,size=sample_size)
index = np.random.uniform(*index_rng,size=sample_size)

sample = np.sort(np.random.randint(0,len(field_data_raw),sample_size))


# Tables to save data into
obj_tbl = Table(names=['objnum','sample_point'],dtype=('i4','i4',))

input_pars = Table(names=inmodel.param_names)

best_pars = Table(names=[pn+'_best'for pn in outmodel.param_names])
mean_pars = Table(names=[pn+'_mean'for pn in outmodel.param_names])
median_pars = Table(names=[pn+'_median'for pn in outmodel.param_names])
mode_pars = Table(names=[pn+'_mode'for pn in outmodel.param_names])


upper_error = Table(names=[pn+'_err_up'for pn in outmodel.param_names])
lower_error = Table(names=[pn+'_err_lo'for pn in outmodel.param_names])
sigma_error = Table(names=[pn+'_stdev'for pn in outmodel.param_names])

# correlation = Table(names=['correlation_matrix'])
lsq = Table(names=['lsq_best'])
rdata = Table(names=['r_data'])

i=0
for s in sample:
	print i,datetime.datetime.now()
	if (i % break_size) == 0 and i>1:
		input_data = hstack([obj_tbl,input_pars])
		output_data = hstack([obj_tbl,best_pars,mean_pars,median_pars,mode_pars,
						sigma_error,upper_error,lower_error,lsq])#,correlation])

		input_data.write(path.join(outdir,"input_data_{:03d}.txt".format(i-break_size)),
						format='ascii.fixed_width')
		output_data.write(path.join(outdir,"output_data_{:03d}.txt".format(i-break_size)),
						format='ascii.fixed_width')


	inmodel.c1_0 = 0.
	inmodel.c2_0 = 0.

	inmodel.E1_2 = E[i,0]
	inmodel.E2_2 = E[i,1]
	inmodel.alpha_2 = A[i]
	inmodel.logI_2 = logI[i]
# 	inmodel.index_2 = index[i]
	inmodel.invindex_2 = index[i]
	
	inmodel.g1_0 = field_data_raw['g1'][s]
	inmodel.g2_0 = field_data_raw['g2'][s]
	inmodel.F1_0 = field_data_raw['F1'][s]
	inmodel.F2_0 = field_data_raw['F2'][s]
	inmodel.G1_0 = field_data_raw['G1'][s]
	inmodel.G2_0 = field_data_raw['G2'][s]
	
	true = inmodel(x,y,p)
	data = true + aimflex.window_image(np.random.randn(*(true.shape)))
	weights = aimflex.window_image(np.ones(true.shape)) # Noise RMS for each pixel...

	### Select only the main image from the generated data image and window down.
	###
	objects = sep.extract(data, 1.5, err=weights)
	
	if len(objects) > 0:
	
		main = np.argmax(objects['flux']) # Pick the bright one
		r=(cut_factor*objects['a'][main]).astype(int)
		xc,yc=objects['x'][main].astype(int), objects['y'][main].astype(int)
		if xc-r < 0:
			r=np.copy(xc)
		if xc+r > data.shape[0]:
			r=data.shape[0] - xc
		if yc-r < 0:
			r=np.copy(yc)
		if yc+r > data.shape[1]:
			r=data.shape[1] - yc
	
		print xc,yc," --- ",r
		print data.shape
		data=aimflex.cut_stamp(data,xc,yc,r)
		true=aimflex.cut_stamp(true,xc,yc,r)
		weights=aimflex.cut_stamp(weights,xc,yc,r)
		print data.shape
	
		objdims = data.shape
		xobj,yobj =np.meshgrid(np.linspace(-(objdims[0]-1)/2,(objdims[0]-1)/2,objdims[0]),
							   np.linspace(-(objdims[1]-1)/2,(objdims[1]-1)/2,objdims[1]))

		# log the offset in the center position
		inmodel.c1_0=xc-dims[0]
		inmodel.c2_0=yc-dims[1]
		true_params = np.copy(inmodel.parameters)

		### FITTING
		###
		aimflex.set_limits(data,outmodel)
		samples, chisqs = aimflex.fit_image(outmodel,data,weights,p,verbose=False)
		samples2D=samples.reshape((-1, outmodel.parameters.size))

		# corner figure
		cfig = corner.corner(samples2D, bins=100, labels=outmodel.param_names,
						truths=np.compress(parmask,true_params), quantiles=[0.16,0.5,0.84])
		cfig.savefig(path.join(outdir,"obj_{:03d}_corner.png".format(i)))
		cfig.clf()
	
		# results to save
		stats = aimflex.sample_stats(samples2D)
		best = np.unravel_index(np.argmin(chisqs),chisqs.shape)
	
		fit_pars = samples[best]
		chisq = chisqs[best]
	
		outmodel.parameters = fit_pars
	
		outimg = outmodel(xobj,yobj,p)

		# Image plots
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
	
		# Parameter track plots
		fig, axes = pl.subplots(outmodel.parameters.size, 1, sharex=True, figsize=(8, 20))

		for k in range(outmodel.parameters.size):
			axes[k].plot(samples2D[:, k].T, color="k", alpha=0.4)
			axes[k].yaxis.set_major_locator(MaxNLocator(5))
			axes[k].axhline(np.compress(parmask,true_params)[k], color="b", lw=2)
			axes[k].set_ylabel(outmodel.param_names[k])

		fig.tight_layout(h_pad=0.0)
		fig.savefig(path.join(outdir,"obj_{:03d}_pars.png".format(i)))
		fig.clf()
	
		del fig, cfig, f

		# Output
		obj_tbl.add_row([i,s])
	
		input_pars.add_row(true_params)
		best_pars.add_row(fit_pars)
		mean_pars.add_row(stats[0])
		median_pars.add_row(stats[2][1])
		mode_pars.add_row(stats[3])
	
		sigma_error.add_row(stats[1])
		upper_error.add_row(stats[2][2] - stats[2][1])
		lower_error.add_row(stats[2][2] - stats[2][1])
	
	# 	correlation.add_row([stats[3]])

		lsq.add_row([chisq])
		rdata.add_row([r])


		pl.close("all")
		i+=1

# Final results
input_data = hstack([obj_tbl,input_pars])
output_data = hstack([obj_tbl,best_pars,mean_pars,median_pars,mode_pars,
					sigma_error,upper_error,lower_error,lsq])#,correlation])

input_data.write(path.join(outdir,"input_data_master.txt"),format='ascii.fixed_width')
output_data.write(path.join(outdir,"output_data_master.txt"),format='ascii.fixed_width')
print 'Finished!!',datetime.datetime.now()

