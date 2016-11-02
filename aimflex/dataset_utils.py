from astropy.io import ascii, fits
import numpy as np

from .core import *


def triptych(data,fit,resid=None,tag='triptych'):
	import matplotlib.pyplot as plt
	"""
		Make a triptych figure with the data, fit, and residual images.
	"""
	if resid is None:
		resid = data - fit
		
	f, (ax1, ax2, ax3) = plt.subplots(1,3)

	ax1.imshow(data,origin='lower')
	ax2.imshow(fit,origin='lower')
	ax3.imshow(resid,origin='lower')

	[a.set_xticks([]) for a in [ax1,ax2,ax3] ]
	[a.set_yticks([]) for a in [ax1,ax2,ax3] ]

	f.savefig(tag+'.png',bbox_inches='tight', transparent=True)
	plt.close()
	


def fit_dataset(image, weights, catalog, outfile, rscale=2.,psf=None,
				ntag='NUMBER',xtag='X_IMAGE',ytag='Y_IMAGE',atag='A_IMAGE',
				outdir=None,save_fig=False, maxiter=10000, maxfun=25000,
				simplex=True):
	"""
		Using a data image, a weight image, and a Source Extractor-like catalog
		fit an AIM profile to each of the objects in the catalog, and store the 
		output.
		
		Assumes the SEcatalog has at least the NUMBER, X_IMAGE, Y_IMAGE, and 
		A_IMAGE columns.
	"""

	sextable = ascii.read(catalog)
	img_hdu = fits.open(image)
	wht_hdu = fits.open(weights)
	
	catno=sextable[ntag]
	cutradii = (sextable[atag]*rscale).astype(int)
	xpix = sextable[xtag].astype(int)
	ypix = sextable[ytag].astype(int)
	
	if outdir is None:
		outdir=""
	elif (outdir[-1] is not "/") and len(outdir)>1:
		outdir = outdir +"/"
	
	
	fit_rows = []
	if simplex:
		print "o hello simplex"
		fitter=AIMSimplexLSQFitter()
	else:
		fitter=AIMLevMarLSQFitter()
	
	for i in range(cutradii.size/40):
		print ("Object %i of %i" % (catno[i],cutradii.size))
		x,y =np.meshgrid(np.linspace(-cutradii[i],cutradii[i],2*cutradii[i]+1),
						 np.linspace(-cutradii[i],cutradii[i],2*cutradii[i]+1))

# 		stamp = window_image(img_hdu[0].data[(ypix[i]-cutradii[i]):(ypix[i]+cutradii[i]+1),\
# 									(xpix[i]-cutradii[i]):(xpix[i]+cutradii[i]+1)])
# 									
# 		weights = window_image(wht_hdu[0].data[(ypix[i]-cutradii[i]):(ypix[i]+cutradii[i]+1),\
# 									 (xpix[i]-cutradii[i]):(xpix[i]+cutradii[i]+1)])
		stamp = cut_stamp(img_hdu[0].data,xpix[i],ypix[i],cutradii[i])
		weights = cut_stamp(wht_hdu[0].data,xpix[i],ypix[i],cutradii[i])

		weights /= np.sum(weights)
		
		outname='obj_%04i_output' % catno[i]
		
		fits.PrimaryHDU(stamp).writeto(outdir+outname+'_stamp.fits',clobber=True)
		fits.PrimaryHDU(weights).writeto(outdir+outname+'_weight.fits',clobber=True)
		
		model=AIM()
# 		model=AIM2()
		set_gaussian_pars(stamp,model,weights=weights)
		initial_cond = model.parameters
		
		if type(model) == AIM:
			model.alpha.max=0.75*cutradii[i]
		if type(model) == AIM2:
			model.M1inv.max = (0.75*cutradii[i])**2
			model.M2inv.max = (0.75*cutradii[i])**2
		model.c1.max=0.5*cutradii[i]
		model.c2.max=0.5*cutradii[i]
		model.c1.min=-0.5*cutradii[i]
		model.c2.min=-0.5*cutradii[i]

		if simplex:
			fit = fitter(model,x,y,psf,stamp,maxiter=maxiter,weights=weights,maxfun=maxfun)
		else:
			# Run MCMC with emcee
			import emcee
			
		print fitter.fit_info.keys()
		print "Initial: ",initial_cond
		print "Final  : ",fit.parameters
		
		fits.PrimaryHDU(fit(x,y,psf)).writeto(outdir+outname+'_fit.fits',clobber=True)
		fits.PrimaryHDU(stamp - fit(x,y,psf)).writeto(outdir+outname+'_residual.fits',clobber=True)
		
		if save_fig:
			triptych(stamp,fit(x,y,psf),resid=(stamp - fit(x,y,psf)),tag=outdir+outname)
		
		outrow=np.append(fit.parameters,
						 [fitter.fit_info['final_func_val'],
						  stamp.size,
						  fitter.fit_info['numiter'],
						  fitter.fit_info['exit_mode']])
		outrow=np.append(outrow,initial_cond)
		fit_rows.append(outrow)
	outcol_names = AIM.param_names+('chisq','npix','niter','exit_mode',)+tuple([x+'_init' for x in AIM.param_names])
	tabledata = Table(rows=fit_rows,names=outcol_names)
	tabledata.write(outdir+outfile,format='ascii')



