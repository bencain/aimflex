import numpy as np
import astropy.modeling as am
import emcee

from astropy.convolution import convolve_fft as cfft
from astropy.modeling.mappings import Identity
from astropy.modeling.fitting import Simplex

from .utilities import *

################# CONSTANTS #########################
N_WALKERS = 100		# emcee walkers
N_BURN = 100			# emcee burn in steps
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
	alpha = am.Parameter(default=5.,min=0.1)
	index = am.Parameter(default=0.5,min=0.01,max=20.)
	E1 = 	am.Parameter(default=0.,min=-E_LIMIT,max=E_LIMIT)
	E2 = 	am.Parameter(default=0.,min=-E_LIMIT,max=E_LIMIT)
	
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


def AIM_lnprior(pars,model):
	"""
		log of the probability 
	"""
	for pn in model.param_names:
		min = getattr(model,pn).min
		max = getattr(model,pn).max
		val = getattr(model,pn).value
		
		if min is not None:
			if val < min:
				return -np.inf
		if max is not None:
			if val > max:
				return -np.inf
	if hasattr(model,'E1_2') and hasattr(model,'E1_2'):
		# IF THE ELLIPTICITY OF THE MODEL PROFILE IS IN
		# THE SHEAR-DEGENERATE PARAMETRIZATION, WE CHECK
		# THE ELLIPTICITY MAGNITUDE FOR COMPLIANCE
		if (model.E1_2**2 + model.E2_2**2) >= E_LIMIT:
			return -np.inf
	return 0.

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
	model.parameters = pars
	lp = AIM_lnprior(pars,model)
	print lp,pars
	if not np.isfinite(lp):
		return -np.inf
	
	return -leastsquare(data, model, weights, x, y, psf)

def fit_image(model,data,weights,psf,verbose=False):
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
		
	p0 = np.repeat(model.parameters,N_WALKERS).reshape((model.parameters.size,N_WALKERS)).transpose()
	sampler = emcee.EnsembleSampler(N_WALKERS, model.parameters.size, 
									AIM_lnprob, 
									args=[model, x, y, psf, data, weights])
	
	# Do the burn in
	pos, prob, state = sampler.run_mcmc(p0, N_BURN)
	sampler.reset()
	
	# Now the main run
	sampler.run_mcmc(pos, N_CHAIN)

	if verbose:
		print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))
		import matplotlib.pyplot as pl
		for i,pn in enumerate(model.param_names):
			pl.ion()
			pl.figure()
			pl.hist(sampler.flatchain[:,i], 100, color="k", histtype="step")
			pl.title("Dimension {}".format(pn))
			pl.ioff()
	
	return sampler.chain
