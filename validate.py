import numpy as np
import aimflex
import matplotlib.pyplot as plt
from astropy.modeling.mappings import Identity
from astropy.io import fits


"""
	This is a validation script to ensure that what goes in comes out.  We'll use the 
	Sersic profile from AIM.
"""

profile=aimflex.sersic
distortion=aimflex.lens_equation
psf = None
ntrial=5
outdir='val_out/'

base_img_axis = 150

pix_scale=0.03 # arcsec/pixel, a la HST

g_lim = 0.5
F_lim = 0.1*pix_scale
G_lim = F_lim
E_lim = 0.7

logI_lim = np.array([-3,0])
alpha_lim = np.array([0,0.5])/pix_scale
index_lim = np.array([0,5.])


g1_in=np.array([])
g2_in=np.array([])
F1_in=np.array([])
F2_in=np.array([])
G1_in=np.array([])
G2_in=np.array([])
E1_in=np.array([])
E2_in=np.array([])

# get input values
while min([g1_in.size,
		   F1_in.size,
		   G1_in.size,
		   E1_in.size]) < ntrial:
	
	if g1_in.size < ntrial:
		g=(np.random.random(2)*2 - 1)*g_lim
		if np.sqrt(np.sum(g**2)) < g_lim:
			g1_in=np.append(g1_in,g[0])
			g2_in=np.append(g2_in,g[1])

	if F1_in.size < ntrial:
		F=(np.random.random(2)*2 - 1)*F_lim
		if np.sqrt(np.sum(F**2)) < F_lim:
			F1_in=np.append(F1_in,F[0])
			F2_in=np.append(F2_in,F[1])

	if G1_in.size < ntrial:
		G=(np.random.random(2)*2 - 1)*G_lim
		if np.sqrt(np.sum(G**2)) < G_lim:
			G1_in=np.append(G1_in,G[0])
			G2_in=np.append(G2_in,G[1])

	if E1_in.size < ntrial:
		E=(np.random.random(2)*2 - 1)*E_lim
		if np.sqrt(np.sum(E**2)) < E_lim:
			E1_in=np.append(E1_in,E[0])
			E2_in=np.append(E2_in,E[1])

logI_in = np.random.random(ntrial)*(logI_lim[1] - logI_lim[0]) + logI_lim[0]
alpha_in = np.random.random(ntrial)*(alpha_lim[1] - alpha_lim[0]) + alpha_lim[0]
index_in = np.random.random(ntrial)*(index_lim[1] - index_lim[0]) + index_lim[0]


# Make images and fit them!
input_pars=[]
output_pars=[]
output_stats=[]

for i in range(ntrial):
	
	input_model = model = (distortion() & Identity(1)) | profile()
	
	input_model.g1_0 = g1_in[i]
	input_model.g2_0 = g2_in[i]
	input_model.F1_0 = F1_in[i]
	input_model.F2_0 = F2_in[i]
	input_model.G1_0 = G1_in[i]
	input_model.G2_0 = G2_in[i]

	input_model.logI_2 = logI_in[i]
	input_model.alpha_2 = alpha_in[i]
	input_model.index_2 = index_in[i]
	input_model.E1_2 = E1_in[i]
	input_model.E2_2 = E2_in[i]

	x,y =np.meshgrid(np.linspace(-base_img_axis,base_img_axis,2*base_img_axis+1),
					 np.linspace(-base_img_axis,base_img_axis,2*base_img_axis+1))

	base_img = input_model(x,y,psf)
	
	gpars = aimflex.parameter_utils.get_gaussian_pars(base_img)
	print gpars
	
	outname='obj_%03i' % i
	print input_model.parameters
	fits.PrimaryHDU(input_model(x,y,psf)).writeto(outdir+outname+'_input.fits',clobber=True)
	
	
	
	
