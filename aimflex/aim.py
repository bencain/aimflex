

import numpy as np
import matplotlib.pyplot as plt

import astropy.modeling as am
from astropy.io import ascii, fits
from astropy.convolution import convolve_fft as cfft
from astropy.table import Table

from astropy.modeling.optimizers import Simplex
# from astropy.modeling.statistic import leastsquare


class AIM(am.FittableModel):
	"""
		This class implements an elliptical Sersic source plane model profile 
		lensed by shear and flexion.  Parameters included:
			logI	- peak surface brightness
			alpha	- size = sqrt(ab) = a*sqrt(E)
			index	- Sersic index (0.5 = Gaussian)
			c1		- image plane x location for beta=0
			c2		- image plane y location for beta=0
			E1		- + polarized ellipticity
			E2		- x polarized ellipticity
			g1		- + polarized reduced shear
			g2		- x polarized reduced shear
			F1		- reduced 1-Flexion x
			F2		- reduced 1-Flexion y
			G1		- reduced 3-Flexion x
			G2		- reduced 3-Flexion y
		We also implement a PSF convolution as well.
	"""
	inputs = ('x','y','psf',)
	outputs = ('img',)
	
	logI  = am.Parameter(default=1.,max=4.)
	alpha = am.Parameter(default=1.,min=0.)
	index = am.Parameter(default=0.5,min=0.01)
	c1 = 	am.Parameter(default=0.)
	c2 = 	am.Parameter(default=0.)
	E1 = 	am.Parameter(default=0.,min=-0.9,max=0.9)
	E2 = 	am.Parameter(default=0.,min=-0.9,max=0.9)
	g1 = 	am.Parameter(default=0.)
	g2 = 	am.Parameter(default=0.)
	F1 = 	am.Parameter(default=0.)
	F2 = 	am.Parameter(default=0.)
	G1 = 	am.Parameter(default=0.)
	G2 = 	am.Parameter(default=0.)
	
	standard_broadcasting = False
		
	@staticmethod
	def evaluate(x,y,psf,
				 logI,alpha,index,c1,c2,
				 E1,E2,g1,g2,F1,F2,G1,G2):	
				
		# Need to get Gaussian2D parameters from my parameters
		if logI > 3.:
			amp=1000.
		else:
			amp = np.power(10.,logI)	# Gaussian amplitude
		
		# Ellipse parameters
		ellipse = convert_epars([alpha,E1,E2],pol_to_ae=True)
		a = ellipse[0]/np.sqrt(ellipse[1])
		b = ellipse[0]*np.sqrt(ellipse[1])
		pa = ellipse[2]
		
		gmodel=am.models.Gaussian2D(amplitude=amp,x_mean=0.,y_mean=0.,
								  x_stddev=a,y_stddev=b,theta=pa)
	
	
		beta0=c1+1j*c2
		E=E1+1j*E2
		g=g1+1j*g2
		F=F1+1j*F2
		G=G1+1j*G2
		
		coo = x+1j*y - beta0
		coo_c = np.conj(coo)
		
		beta = coo - g*coo_c - np.conj(F)*coo**2 \
				- 2*F*coo*coo_c - G*coo_c**2 
				
		# Now for the model:
		img = gmodel(beta.real,beta.imag)
		out = img * 0.0
		
		if np.any(np.isnan(psf)):
			out = img
		elif psf is None:
			out = img
		else:
			if len(psf.shape) == 2:
				out = cfft(img,psf)
		
		## WINDOW THE IMAGE - note on assumptions: x,y are zero centered & have odd width/height
		out[(2*x/(x.shape[0]-1))**2+(2*y/(x.shape[1]-1))**2>1] = 0.
		return out

class AIM2(am.FittableModel):
	"""
		This class implements an elliptical Sersic source plane model profile 
		lensed by shear and flexion.  Parameters included:
			logI	- peak surface brightness
			index	- Sersic index (0.5 = Gaussian)
			c1		- image plane x location for beta=0
			c2		- image plane y location for beta=0
			1_M1	- coefficient of the x**2 term in the Gaussian argument
			1_M2	- coefficient of the y**2 term in the Gaussian argument
			M3		- coefficient of the x*y term in the Gaussian argument
			g1		- + polarized reduced shear
			g2		- x polarized reduced shear
			F1		- reduced 1-Flexion x
			F2		- reduced 1-Flexion y
			G1		- reduced 3-Flexion x
			G2		- reduced 3-Flexion y
		We also implement a PSF convolution as well.
	"""
	inputs = ('x','y','psf',)
	outputs = ('img',)
	
	logI  = am.Parameter(default=1.,max=4.)
	index = am.Parameter(default=0.5,min=0.01)
	c1 = 	am.Parameter(default=0.)
	c2 = 	am.Parameter(default=0.)
	M1inv = am.Parameter(default=5.,min=1e-3)
	M2inv = am.Parameter(default=5.,min=1e-3)
	M3 = 	am.Parameter(default=0.)
	g1 = 	am.Parameter(default=0.)
	g2 = 	am.Parameter(default=0.)
	F1 = 	am.Parameter(default=0.)
	F2 = 	am.Parameter(default=0.)
	G1 = 	am.Parameter(default=0.)
	G2 = 	am.Parameter(default=0.)
	
	standard_broadcasting = False
		
	@staticmethod
	def evaluate(x,y,psf,
				 logI,index,c1,c2,
				 M1inv,M2inv,M3,g1,g2,F1,F2,G1,G2):	
				
		aee = convert_epars([1/M1inv, 1/M2inv, M3],mi_to_pol=True)
		
		return AIM.evaluate(x,y,psf,logI,aee[0],index,c1,c2,aee[1],aee[2],
							g1,g2,F1,F2,G1,G2)


def _convert_input(x, y, p, z=None, n_models=1, model_set_axis=0):
	"""Convert inputs to float arrays."""

	x = np.asarray(x, dtype=np.float)
	y = np.asarray(y, dtype=np.float)
	p = np.asarray(p, dtype=np.float)
	if z is not None:
		z = np.asarray(z, dtype=np.float)

	# For compatibility with how the linear fitter code currently expects to
	# work, shift the dependent variable's axes to the expected locations
	if n_models > 1:
		if z is None:
			if y.shape[model_set_axis] != n_models:
				raise ValueError(
					"Number of data sets (y array is expected to equal "
					"the number of parameter sets)")
			# For a 1-D model the y coordinate's model-set-axis is expected to
			# be last, so that its first dimension is the same length as the x
			# coordinates.  This is in line with the expectations of
			# numpy.linalg.lstsq:
			# http://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.lstsq.html
			# That is, each model should be represented by a column.  TODO:
			# Obviously this is a detail of np.linalg.lstsq and should be
			# handled specifically by any fitters that use it...
			y = np.rollaxis(y, model_set_axis, y.ndim)
		else:
			# Shape of z excluding model_set_axis
			z_shape = z.shape[:model_set_axis] + z.shape[model_set_axis + 1:]

			if not (x.shape == y.shape == z_shape):
				raise ValueError("x, y and z should have the same shape")

	if z is None:
		farg = (x, y, p)
	else:
		farg = (x, y, p, z)
	
	return farg

def leastsquare(measured_vals, updated_model, weights, x, y=None, p=None):
    """
    Least square statistic with optional weights.

    Parameters
    ----------
    measured_vals : `~numpy.ndarray`
        Measured data values.
    updated_model : `~astropy.modeling.Model`
        Model with parameters set by the current iteration of the optimizer.
    weights : `~numpy.ndarray`
        Array of weights to apply to each residual.
    x : `~numpy.ndarray`
        Independent variable "x" to evaluate the model on.
    y : `~numpy.ndarray`, optional
        Independent variable "y" to evaluate the model on, for 2D models.
    p : `~numpy.ndarray`, optional
    	Fixed model parameters, such as a psf for a 2D image model

    Returns
    -------
    res : float
        The sum of least squares.
    """

    if y is None and p is None:
        model_vals = updated_model(x)
    else:
        model_vals = updated_model(x, y, p)
    if weights is None:
        return np.sum((model_vals - measured_vals) ** 2)
    else:
        return np.sum((weights * (model_vals - measured_vals)) ** 2)

class AIMSimplexLSQFitter(am.fitting.SimplexLSQFitter):
    """

    Simplex algorithm and least squares statistic.

    Raises
    ------
    ModelLinearityError
        A linear model is passed to a nonlinear fitter

    """

    supported_constraints = Simplex.supported_constraints

    def __init__(self):
        super(am.fitting.SimplexLSQFitter, self).__init__(optimizer=Simplex,
                                                  statistic=leastsquare)
        self.fit_info = {}



    def __call__(self, model, x, y, p, z=None, weights=None, **kwargs):
        """
        Fit data to this model.

        Parameters
        ----------
        model : `~astropy.modeling.FittableModel`
            model to fit to x, y, p, z
        x : array
            input coordinates
        y : array
            input coordinates
        p : array
        	additional parameter input (e.g., PSF)
        z : array (optional)
            input coordinates
        weights : array (optional)
            weights
        kwargs : dict
            optional keyword arguments to be passed to the optimizer or the statistic  		

        maxiter : int
            maximum number of iterations
        epsilon : float
            the step size for finite-difference derivative estimates
        acc : float
            Relative error in approximate solution

        Returns
        -------
        model_copy : `~astropy.modeling.FittableModel`
            a copy of the input model with parameters set by the fitter
        """

        model_copy = am.fitting._validate_model(model,
                                     self._opt_method.supported_constraints)
        
        farg = _convert_input(x, y, p, z,n_models=len(model_copy))
        farg = (model_copy, weights, ) + farg

        p0, _ = am.fitting._model_to_fit_params(model_copy)
		
        fitparams, self.fit_info = self._opt_method(self.objective_function, p0, farg, **kwargs)
        
        am.fitting._fitter_to_model_params(model_copy, fitparams)
        return model_copy




def set_gaussian_pars(image,model,weights=1.):
	"""
		A function for creating Gaussian model initial parameters from the data and
		setting them into the model object.
	"""
	
	cts = calc_moments(image,weights=weights)
	ctr = [calc_moments(image,xord=1,weights=weights)/cts,
		   calc_moments(image,yord=1,weights=weights)/cts]
	Q11 = calc_moments(image,xord=2,weights=weights)/cts
	Q22 = calc_moments(image,yord=2,weights=weights)/cts
	Q12 = calc_moments(image,xord=1,yord=1,weights=weights)/cts
	
	m1=Q11-Q12**2/Q22
	m2=Q22-Q12**2/Q11
	m3=-Q12/(Q11*Q22-Q12**2)
	
	epars = convert_epars([m1,m2,m3],mi_to_pol=True)
	
	logI=np.log10(cts/(2*np.pi*epars[0]**2))
	
	model.logI = logI
	model.c1=ctr[0]
	model.c2=ctr[1]
	model.alpha=epars[0]
	if hasattr(model,'E1'):
		model.E1=epars[1]
		model.E2=epars[2]
	elif hasattr(model,'S1'):
		model.S1=epars[1]
		model.S2=epars[2]

	
	
	
def STUV(E,g,F,G):
	"""
		Convert ellipticity and lensing parameters to the combined 
		lensing parametrization.
	"""
	M = 1 + E*np.conj(g)
	S = (E+g)/M
	T = (F - np.conj(E)*G)/np.conj(M)
	U = (F - E*np.conj(F))/M
	V = (G - E*F)/M
	
	return S,T,U,V



def fit_dataset(image, weights, catalog, outfile, rscale=2.,psf=None,
				ntag='NUMBER',xtag='X_IMAGE',ytag='Y_IMAGE',atag='A_IMAGE',
				outdir=None,save_fig=False, maxiter=10000, maxfun=25000):
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
	fitter=AIMSimplexLSQFitter()
	
	for i in range(cutradii.size/40):
		print ("Object %i of %i" % (catno[i],cutradii.size))
		x,y =np.meshgrid(np.linspace(-cutradii[i],cutradii[i],2*cutradii[i]+1),
						 np.linspace(-cutradii[i],cutradii[i],2*cutradii[i]+1))
		mask=(x**2+y**2<=cutradii[i]**2).astype(int)
		
		stamp = img_hdu[0].data[(ypix[i]-cutradii[i]):(ypix[i]+cutradii[i]+1),\
									(xpix[i]-cutradii[i]):(xpix[i]+cutradii[i]+1)]*mask
									
		weights = wht_hdu[0].data[(ypix[i]-cutradii[i]):(ypix[i]+cutradii[i]+1),\
									 (xpix[i]-cutradii[i]):(xpix[i]+cutradii[i]+1)]*mask
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

		fit = fitter(model,x,y,psf,stamp,maxiter=maxiter,weights=weights,maxfun=maxfun)
		print fitter.fit_info
				
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
	
#####################################################################################
#################	UTILITIES	#####################################################
#####################################################################################
	
def calc_moments(image,xord=0,yord=0,weights=1.):
	"""
		A function for calculating arbitrary xy image moments
	"""
	dims=image.shape
	
	# Handle weights
	if not type(weights) == np.ndarray:
		w = 1.
	else:
		w = weights/np.sum(weights) # Normalize

	if xord < 0: # No negative orders
		xord=0
	if yord < 0:
		yord=0
	
	x,y =np.meshgrid(np.linspace(-(dims[0]-1)/2,(dims[0]-1)/2,dims[0]),
					 np.linspace(-(dims[1]-1)/2,(dims[1]-1)/2,dims[1]))
	
	if (xord+yord) > 1: # Calc ctr and flux for higher moments
		cts = calc_moments(image,xord=0,yord=0,weights=w)
		ctr = [calc_moments(image,xord=1,yord=0,weights=w)/cts,
			   calc_moments(image,xord=0,yord=1,weights=w)/cts]
	else:
		ctr=[0.,0.]
	
	return np.sum(image*w*np.power((x-ctr[0]),xord)*np.power((y-ctr[1]),yord))
	
	

def convert_epars(epars,
                  ae_to_pol=False, pol_to_mi=False,
                  mi_to_ae=False,  pol_to_ae=False,
                  mi_to_pol=False, ae_to_mi=False):
	"""
		This function converts between parametrizations of ellipse
		parameters. We consider an elliptical profile given by
			I(x,y)=I(r)
			
			r^2 = ((x-X_c)*cos(xi)+(y-Y_c)*sin(xi))^2/A^2 +
				  ((y-Y_c)*cos(xi)-(x-X_c)*sin(xi))^2/(eps*A)^2

		The center (X_c,Y_c) is not included in this function, just the 
		three necessary shape parameters defined by one of the following 
		parametrizations:

		AE: Semimajor/minor axis geometric mean (AB)^0.5,
			axis ratio EPS and position angle XI

		POL: alpha=(AB)^0.5, "plus-mode" ellipticity EPLUS and
			 "cross-mode" ellipticity ECROSS (evoking polarization)
			 
			 r^2=(1/alpha^2) * (1/sqrt(1-eplus^2-ecross^2) *
				((1-eplus)*(x-X_c)^2 + (1+eplus)*(y-Y_c) -
				2*ecross*(x-X_c)*(y-Y_c))

		MI:  Inverse x^2 coefficient M1, inverse y^2 coefficient M2 and xy
			 coefficient M3 from
			 r^2 = ((x-X_c)^2/M1^2 + (y-Y_c)^2/M2 + 2*M3*(x-X_c)*(y-Y_c))
	"""
	outpars=np.copy(epars)
	if np.isscalar(epars):
		print "aim.convert_epars requires a 3 element iterable input."
		print " -> this looks like a scalar: ",epars
		return outpars
	
	if len(epars) < 3:
		print "aim.convert_epars requires a 3 element iterable input."
		print " -> this looks like too few: ",epars
		return outpars
	elif len(epars) > 3:
		print "aim.convert_epars requires a 3 element iterable input."
		print " -> this looks like too many: ",epars
		return outpars
	
	conds=np.array([ae_to_pol, pol_to_mi, mi_to_ae, 
					pol_to_ae, mi_to_pol, ae_to_mi],dtype=bool)
					
	if len(conds[conds]) == 0:
		print "Doesn't look like a conversion was selected."
		return outpars
	elif len(conds[conds]) > 1:
		print "Looks like multiple conversions were selected."
		return outpars


	if ae_to_pol:
		A=epars[0]
		eps=epars[1]
		xi=epars[2]
		
		alpha=A*np.sqrt(eps)
		
		# This is
		# E = ( Q11 - Q22 + 2iQ12 )/( Q11 + Q22 )
		# eplus= ((1d)-eps^2)*cos((2d)*xi)/((1d)+eps^2)
		# ecross=((1d)-eps^2)*sin((2d)*xi)/((1d)+eps^2)
		
		# E = ( Q11 - Q22 + 2iQ12 )/( Q11 + Q22 +2sqrt( Q11*Q22 - Q12^2 ) )
		
		eplus= (1.-eps)*np.cos(2.*xi)/(1.+eps)
		ecross=(1.-eps)*np.sin(2.*xi)/(1.+eps)
		
		outpars[0:3]=alpha,eplus,ecross

	elif pol_to_mi:
		alpha=epars[0]
		eplus=epars[1]
		ecross=epars[2]

		emag=np.sqrt(eplus**2+ecross**2)

		# This is
		# E = ( Q11 - Q22 + 2iQ12 )/( Q11 + Q22 )
		# m1=alpha^2*sqrt(1d -emag^2)/(1d -eplus)
		# m2=alpha^2*sqrt(1d -emag^2)/(1d +eplus)
		# m3=-ecross/(alpha^2*sqrt(1d -emag^2))

		# E = ( Q11 - Q22 + 2iQ12 )/( Q11 + Q22 +2sqrt( Q11*Q22 - Q12^2 ) )
		m1=alpha**2*(1. - emag**2)/(1. + emag**2 - 2.*eplus)
		m2=alpha**2*(1. - emag**2)/(1. + emag**2 + 2.*eplus)
		m3=-2.*ecross/(alpha**2*(1. -emag**2))

		outpars[0:3]=m1,m2,m3

	elif pol_to_ae:
		alpha=epars[0]
		eplus=epars[1]
		ecross=epars[2]
		emag=np.sqrt(eplus**2+ecross**2)

		# This is
		# E = ( Q11 - Q22 + 2iQ12 )/( Q11 + Q22 )
		# eps=sqrt((1d -emag)/(1d +emag))
		# A=alpha/sqrt(eps)
		# xi=0.5d*atan(ecross,eplus)
		
		# E = ( Q11 - Q22 + 2iQ12 )/( Q11 + Q22 +2sqrt( Q11*Q22 - Q12^2 ) )
		eps=(1. - emag)/(1. + emag)
		A=alpha/np.sqrt(eps)
		xi=0.5*np.arctan2(ecross,eplus)

		outpars[0:3]=A,eps,xi

	elif mi_to_pol:
		m1=epars[0]
		m2=epars[1]
		m3=epars[2]

		# This is
		# E = ( Q11 - Q22 + 2iQ12 )/( Q11 + Q22 )
		# alpha=(m1*m2/(1d - m1*m2*m3^2))^(0.25d)
		# eplus=(m1-m2)/(m1+m2)
		# ecross=-2d*m1*m2*m3/(m1+m2)

		# E = ( Q11 - Q22 + 2iQ12 )/( Q11 + Q22 +2sqrt( Q11*Q22 - Q12^2 ) )
		alpha=np.power(m1*m2/(1. - m1*m2*m3**2),0.25)
		eplus=(m1 - m2)/(m1 + m2 + 2.*np.sqrt(m1*m2 - (m1*m2*m3)**2))
		ecross=-2.*m1*m2*m3/(m1 + m2 + 2*np.sqrt(m1*m2 - (m1*m2*m3)**2))

		outpars[0:3]=alpha,eplus,ecross

	elif ae_to_mi:
		outpars=convert_epars(convert_epars(epars,ae_to_pol=True),pol_to_mi=True)

	elif mi_to_ae:
		outpars=convert_epars(convert_epars(epars,mi_to_pol=True),pol_to_ae=True)

	return outpars


def triptych(data,fit,resid=None,tag='triptych'):
	"""
		Make a triptych figure with the data, fit, and residual images.
	"""
	if resid is None:
		resid = data - fit
		
	f, (ax1, ax2, ax3) = plt.subplots(1,3)

	ax1.imshow(data)
	ax2.imshow(fit)
	ax3.imshow(resid)

	[a.set_xticks([]) for a in [ax1,ax2,ax3] ]
	[a.set_yticks([]) for a in [ax1,ax2,ax3] ]

	f.savefig(tag+'.png',bbox_inches='tight', transparent=True)
	plt.close()
	
