import aimflex
import corner
import sep
import datetime
import sys

import numpy as np
import matplotlib.pyplot as pl

from matplotlib.ticker import MaxNLocator
from matplotlib.patches import Ellipse
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

sample_size=10
break_size =10

dims = (201,201)
x,y =np.meshgrid(np.linspace(-(dims[0]-1)/2,(dims[0]-1)/2,dims[0]),
                 np.linspace(-(dims[1]-1)/2,(dims[1]-1)/2,dims[1]))
p=None
model = aimflex.get_aim(gaussian=True)

E = np.random.randn(sample_size,2)*e_sigma
A = np.random.uniform(*alpha_rng,size=sample_size)
logI = np.random.uniform(*logI_rng,size=sample_size)
index = np.random.uniform(*index_rng,size=sample_size)

sample = np.sort(np.random.randint(0,len(field_data_raw),sample_size))

for i,s in enumerate(sample):
	model.c1_0 = 0.
	model.c2_0 = 0.

	model.E1_2 = E[i,0]
	model.E2_2 = E[i,1]
	model.alpha_2 = A[i]
	model.logI_2 = logI[i]
	model.invindex_2 = index[i]
	
	model.g1_0 = field_data_raw['g1'][s]
	model.g2_0 = field_data_raw['g2'][s]
	model.F1_0 = field_data_raw['F1'][s]
	model.F2_0 = field_data_raw['F2'][s]
	model.G1_0 = field_data_raw['G1'][s]
	model.G2_0 = field_data_raw['G2'][s]
	
	true_params = np.copy(model.parameters)
	
	true = model(x,y,p)
	data = true + aimflex.window_image(np.random.randn(*(true.shape)))
	weights = aimflex.window_image(np.ones(true.shape))

	objects = sep.extract(data, 1.5, err=weights)
	
	main = np.argmax(objects['flux'])
	r=(3.0*objects['a'][main]).astype(int)
	xc,yc=objects['x'][main].astype(int), objects['y'][main].astype(int)
	
	print xc,yc,r
	cutdata=aimflex.cut_stamp(data,xc,yc,r)
	
	# plot background-subtracted image
	fig, ax = pl.subplots()
	im = ax.imshow(cutdata, interpolation='nearest',origin='lower')

# plot an ellipse for each object
	for i in range(len(objects)):
		e = Ellipse(xy=(objects['x'][i]-xc+r, objects['y'][i]-yc+r),
					width=6*objects['a'][i],
					height=6*objects['b'][i],
					angle=objects['theta'][i] * 180. / np.pi)
		e.set_facecolor('none')
		e.set_edgecolor('red')
		ax.add_artist(e)
	pl.show()
