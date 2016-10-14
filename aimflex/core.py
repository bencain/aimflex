import numpy as np
import astropy.modeling as am
import emcee

from astropy.convolution import convolve_fft as cfft
from astropy.modeling.mappings import Identity

from .utilities import *

################# CONSTANTS #########################
N_WALKERS = 250		# emcee walkers
N_BURN = 100		# emcee burn in steps
N_CHAIN = 1000		# emcee MCMC steps per walker

E_LIMIT = 0.9		# Ellipticity magnitude limit
################# CLASSES ###########################

class lens_equation(am.FittableModel):
	"""
		This class implements a quadratic local lens distortion with
		reduced shear, and flexion as inputs.  Parameters are
			c1		- image plane x location for beta=0
			c2		- image plane y location for beta=0
			g1		- + polarized reduced shear
			g2		- x polarized reduced shear
			F1		- reduced 1-Flexion x
			F2		- reduced 1-Flexion y
			G1		- reduced 3-Flexion x
			G2		- reduced 3-Flexion y
		We assume a complex image plane position array as input, and
		return the same for the source plane as output.  Flexion units
		are assumed to be inverse of the image/source plane units.
	"""
	inputs = ('theta1','theta2',)
	outputs = ('beta1','beta2',)
	
	c1 = 	am.Parameter(default=0.)
	c2 = 	am.Parameter(default=0.)
	g1 = 	am.Parameter(default=0.)
	g2 = 	am.Parameter(default=0.)
	F1 = 	am.Parameter(default=0.)
	F2 = 	am.Parameter(default=0.)
	G1 = 	am.Parameter(default=0.)
	G2 = 	am.Parameter(default=0.)
	
# 	standard_broadcasting = False
		
	@staticmethod
	def evaluate(theta1,theta2, c1,c2, g1,g2, F1,F2,G1,G2):	
		
		theta=theta1+1j*theta2
		theta0=c1+1j*c2
		
		g=g1+1j*g2
		F=F1+1j*F2
		G=G1+1j*G2
		
		coo = theta - theta0
		coo_c = np.conj(coo)
		
		beta = coo - g*coo_c - np.conj(F)*coo**2 \
				- 2*F*coo*coo_c - G*coo_c**2 				

		return beta.real,beta.imag


class sersic(am.FittableModel):
	"""
		This class implements an elliptical Sersic source plane model profile.
		Parameters included:
			logI	- peak surface brightness
			alpha	- size = sqrt(ab) = a*sqrt(E)
			index	- Sersic index (0.5 = Gaussian)
			E1		- + polarized ellipticity
			E2		- x polarized ellipticity
		We also implement a PSF convolution as well.
	"""
	inputs = ('x','y','psf',)
	outputs = ('img',)

	logI  = am.Parameter(default=1.,max=4.)
	alpha = am.Parameter(default=1.,min=0.)
	index = am.Parameter(default=0.5,min=0.01,max=20.)
	E1 = 	am.Parameter(default=0.,min=-1,max=1)
	E2 = 	am.Parameter(default=0.,min=-1,max=1)
	
	standard_broadcasting = False
		
	@staticmethod
	def evaluate(x,y,psf,
				 logI,alpha,index,E1,E2):	
				
		# Need to get Gaussian2D parameters from my parameters
		if logI > logI.max:
			logI= logI.max
		
		amp = np.power(10.,logI)	# Gaussian amplitude
		
		# Ensure ellipticities don't get out of whack
		emag=np.sqrt(E1**2+E2**2)
		if emag > 1:
			E1/=emag
			E2/=emag
		
		# Ellipse parameters
		ellipse = convert_epars([alpha,E1,E2],pol_to_ae=True)
		a = ellipse[0]
		b = ellipse[0]*ellipse[1]
		pa = ellipse[2]

						
		# Now for the model:
		beta = x+1j*y
		
		beta *= np.exp(-1j*pa)
		r = np.sqrt((beta.real/a)**2 + (beta.imag/b)**2)
		img = amp*np.exp(-np.power(r,1./index))
		out = img * 0.0
		
		if np.any(np.isnan(psf)):
			out = img
		elif psf is None:
			out = img
		else:
			if len(psf.shape) == 2:
				out = cfft(img,psf)
		
		## WINDOW THE IMAGE 
		return window_image(out)



################# METHODS ###########################

def get_AIM():
	"""
		Return a lensed Sersic model object.
	"""
	
	profile=sersic
	distortion=lens_equation
	return (distortion() & Identity(1)) | profile()


def leastsquare(measured_vals, model, weights, x, y=None, p=None):
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

    if y is None:
        model_vals = model(x)
    else:
        model_vals = model(x, y, p)
    if weights is None:
        return np.sum((model_vals - measured_vals) ** 2)
    else:
        return np.sum((weights * (model_vals - measured_vals)) ** 2)


def AIM_lnprob(pars, model, x, y, psf, data, weights):
	"""
		Log of the probability distribution, using the leastsquare 
		function as a chis-squared estimator for the likelihood.
		
		We limit the max and min values based on the settings of the
		model parameters, and if there is an AIM style setup, i.e., a
		lens model distortion passed to a source plane profile, and if
		that profile has E1 and E2 parameters (degenerate with shear)
		then we restrict those parameters to be less than the global 
		E_LIMIT
	"""
	for pn in model.parameter_names:
		min = getattr(model,pn).min
		max = getattr(model,pn).max
		val = getattr(model,pn).value
		
		if min is not None:
			if val < min:
				return -1e9
		if max is not None:
			if val > max:
				return -1e9
	if hasattr(model,'E1_2') and hasattr(model,'E1_2'):
		# IF THE ELLIPTICITY OF THE MODEL PROFILE IS IN
		# THE SHEAR-DEGENERATE PARAMETRIZATION, WE CHECK
		# THE ELLIPTICITY MAGNITUDE FOR COMPLIANCE
		if (model.E1_2**2 + model.E2_2**2) >= E_LIMIT:
			return -1e9

	model.parameters = pars
	return -leastsquare(data, model, weights, x, y, psf)

def fit_image(model,data,weights,psf):
	"""
		Use emcee to find the best parameter fit plus errors. Takes as
		input:
			
			model - an AIM model object compatible with AIM_lnprob
			
			data - a 2D data image with odd pixel dimensions
			
			weights - a 2D weight image of the same shape as data
			
			psf - a 2D point-spread function image of the same pixel
				  scale as data and weights
		
		We assume that the model has already been set to the initial guess
		paramter values.
	"""
	
	# Set up the x and y grids
	axes = [(s-1)/2. for s in data.shape[::-1]]
	x,y =np.meshgrid(np.linspace(-axes[0],axes[0],2*axes[0]+1),
					 np.linspace(-axes[1],axes[1],2*axes[1]+1))
	
	p0 = model.parameters
	sampler = emcee.EnsembleSampler(nwalkers, model.parameters.size, 
									AIM_lnprob, 
									args=[model, x, y, psf, data, weights])
	
	# Do the burn in
	pos, prob, state = sampler.run_mcmc(p0, N_BURN)
	sampler.reset()
	
	# Now the main run
	sampler.run_mcmc(pos, N_CHAIN)
	
	