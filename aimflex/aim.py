import numpy as np
import astropy.modeling as am
from astropy.convolution import convolve_fft as cfft
from astropy.table import Table
from astropy.modeling.optimizers import Simplex
from astropy.modeling.mappings import Identity

from .parameter_utils import *

DEFAULT_MAXITER = 100
DEFAULT_ACC = 1e-07
DEFAULT_EPS = np.sqrt(np.finfo(float).eps)


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
	index = am.Parameter(default=0.5,min=0.01,max=20.)
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
# 		a = ellipse[0]/np.sqrt(ellipse[1]) this is wrong
# 		b = ellipse[0]*np.sqrt(ellipse[1])
		a = ellipse[0]
		b = ellipse[0]*ellipse[1]
		pa = ellipse[2]
		
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



class AIMLevMarLSQFitter(object):
    """
    !!!!!!!!!!!-----------MODIFIED FOR AIM--------------!!!!!!!!!!!
    Based on the astropy implementation (heavily)
    
    Levenberg-Marquardt algorithm and least squares statistic.

    Attributes
    ----------
    fit_info : dict
        The `scipy.optimize.leastsq` result for the most recent fit (see
        notes).

    Notes
    -----
    The ``fit_info`` dictionary contains the values returned by
    `scipy.optimize.leastsq` for the most recent fit, including the values from
    the ``infodict`` dictionary it returns. See the `scipy.optimize.leastsq`
    documentation for details on the meaning of these values. Note that the
    ``x`` return value is *not* included (as it is instead the parameter values
    of the returned model).

    Additionally, one additional element of ``fit_info`` is computed whenever a
    model is fit, with the key 'param_cov'. The corresponding value is the
    covariance matrix of the parameters as a 2D numpy array.  The order of the
    matrix elements matches the order of the parameters in the fitted model
    (i.e., the same order as ``model.param_names``).
    """

    supported_constraints = ['fixed', 'tied', 'bounds']
    """
    The constraint types supported by this fitter type.
    """

    def __init__(self):
        self.fit_info = {'nfev': None,
                         'fvec': None,
                         'fjac': None,
                         'ipvt': None,
                         'qtf': None,
                         'message': None,
                         'ierr': None,
                         'param_jac': None,
                         'param_cov': None,
                         'final_func_val':None}

        super(AIMLevMarLSQFitter, self).__init__()

    def objective_function(self, fps, *args):
        """
        Function to minimize.

        Parameters
        ----------
        fps : list
            parameters returned by the fitter
        args : list
            [model, [weights], [input coordinates]]
        """

        model = args[0]
        weights = args[1]
        am.fitting._fitter_to_model_params(model, fps)
        meas = args[-1]
        if weights is None:
            return np.ravel(model(*args[2 : -1]) - meas)
        else:
            return np.ravel(weights * (model(*args[2 : -1]) - meas))

    def __call__(self, model, x, y, p, z=None, weights=None,
                 maxiter=DEFAULT_MAXITER, acc=DEFAULT_ACC,
                 epsilon=DEFAULT_EPS, estimate_jacobian=False):

        """
        Fit data to this model.

        Parameters
        ----------
        model : `~astropy.modeling.FittableModel`
            model to fit to x, y, z
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
        maxiter : int
            maximum number of iterations
        acc : float
            Relative error desired in the approximate solution
        epsilon : float
            A suitable step length for the forward-difference
            approximation of the Jacobian (if model.fjac=None). If
            epsfcn is less than the machine precision, it is
            assumed that the relative errors in the functions are
            of the order of the machine precision.
        estimate_jacobian : bool
            If False (default) and if the model has a fit_deriv method,
            it will be used. Otherwise the Jacobian will be estimated.
            If True, the Jacobian will be estimated in any case.

        Returns
        -------
        model_copy : `~astropy.modeling.FittableModel`
            a copy of the input model with parameters set by the fitter
        """

        from scipy import optimize

        model_copy = am.fitting._validate_model(model, self.supported_constraints)
        farg = (model_copy, weights, ) + _convert_input(x, y, p, z)

        if model_copy.fit_deriv is None or estimate_jacobian:
            dfunc = None
        else:
            dfunc = self._wrap_deriv
        init_values, _ = am.fitting._model_to_fit_params(model_copy)
        fitparams, cov_x, dinfo, mess, ierr = optimize.leastsq(
            self.objective_function, init_values, args=farg, Dfun=dfunc,
            col_deriv=model_copy.col_fit_deriv, maxfev=maxiter, epsfcn=epsilon,
            xtol=acc, full_output=True)
        am.fitting._fitter_to_model_params(model_copy, fitparams)
        self.fit_info.update(dinfo)
        self.fit_info['cov_x'] = cov_x
        self.fit_info['message'] = mess
        self.fit_info['ierr'] = ierr
        if ierr not in [1, 2, 3, 4]:
            warnings.warn("The fit may be unsuccessful; check "
                          "fit_info['message'] for more information.",
                          AstropyUserWarning)

        # now try to compute the true covariance matrix
        if (len(y) > len(init_values)) and cov_x is not None:
            sum_sqrs = np.sum(self.objective_function(fitparams, *farg)**2)
            dof = len(y) - len(init_values)
            self.fit_info['param_cov'] = cov_x * sum_sqrs / dof
        else:
            self.fit_info['param_cov'] = None

        return model_copy

    @staticmethod
    def _wrap_deriv(params, model, weights, x, y, z=None):
        """
        Wraps the method calculating the Jacobian of the function to account
        for model constraints.

        `scipy.optimize.leastsq` expects the function derivative to have the
        above signature (parlist, (argtuple)). In order to accommodate model
        constraints, instead of using p directly, we set the parameter list in
        this function.
        """

        if weights is None:
            weights = 1.0

        if any(model.fixed.values()) or any(model.tied.values()):

            if z is None:
                full_deriv = np.ravel(weights) * np.array(model.fit_deriv(x, *model.parameters))
            else:
                full_deriv = (np.ravel(weights) * np.array(model.fit_deriv(x, y, *model.parameters)).T).T

            pars = [getattr(model, name) for name in model.param_names]
            fixed = [par.fixed for par in pars]
            tied = [par.tied for par in pars]
            tied = list(np.where([par.tied is not False for par in pars],
                                 True, tied))
            fix_and_tie = np.logical_or(fixed, tied)
            ind = np.logical_not(fix_and_tie)

            if not model.col_fit_deriv:
                full_deriv = np.asarray(full_deriv).T
                residues = np.asarray(full_deriv[np.nonzero(ind)]).T
            else:
                residues = full_deriv[np.nonzero(ind)]

            return [np.ravel(_) for _ in residues]
        else:
            if z is None:
                return [np.ravel(_) for _ in np.ravel(weights) * np.array(model.fit_deriv(x, *params))]
            else:
                return [np.ravel(_) for _ in (np.ravel(weights) * np.array(model.fit_deriv(x, y, *params)).T).T]




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


def window_image(image):
	"""
		This function windows an image such that all pixels outside the 
		ellipse bounded by the edges of the image and whose major and minor 
		axes align with the axes of the image.
	"""
	## REMEMBER - PYTHON STORES INDICES Y,X
	axes = [(s-1)/2. for s in image.shape[::-1]]
	
	# Need some testing to stop bad things from happening if input is bad 
	# if it or doesn't have shape...eventually
	
	x,y =np.meshgrid(np.linspace(-axes[0],axes[0],2*axes[0]+1),
					 np.linspace(-axes[1],axes[1],2*axes[1]+1))
	out = np.copy(image)
	out[(x/axes[0])**2 + (y/axes[1])**2 > 1]=0.
	
	return out
	


	
#####################################################################################
#################	new version	#####################################################
#####################################################################################
	



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
	E1 = 	am.Parameter(default=0.,min=-0.9,max=0.9)
	E2 = 	am.Parameter(default=0.,min=-0.9,max=0.9)
	
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


def finite_diff_hessian(model,data,*argpars,**kwargs):
	
	if 'step' not in kwargs:
		step=1e-6
	elif kwargs['step'] is None:
		step=1e-6
	else:
		step=kwargs['step']
	
	if 'weights' not in kwargs:
		weights=1
	elif kwargs['weights'] is None:
		weights=1
	else:
		weights=kwargs['weights']
	
	model_copy = model.copy()
	pnames=np.copy(model_copy.param_names)
	pvals =np.copy(model_copy.parameters)
	
	h = np.array(step)
	if h.size == 1:
		h = h*np.ones(pvals.size)
	if h.size < pvals.size:
		h = h[0]*np.ones(pvals.size)
	
	def fcn1(h,k):
		model_copy.parameters=pvals
		model_copy.parameters[k]=pvals[k]+h
		return np.sum(weights*(model_copy(*argpars) - data)**2)
	def fcn2(h1,h2,k1,k2):
		model_copy.parameters=pvals
		model_copy.parameters[k1]=pvals[k1]+h1
		model_copy.parameters[k2]=pvals[k2]+h2
		return np.sum(weights*(model_copy(*argpars) - data)**2)
	
	
	hess = np.zeros((pvals.size,pvals.size))
	for i in range(pvals.size):
		for j in np.arange(i,pvals.size):
			if i==j:
				hess[i,j] = (-fcn1(2*h[i],i) + 16*fcn1(h[i],i) - 30*fcn1(0,i) + 16*fcn1(-h[i],i) - fcn1(-2*h[i],i))/(12*h[i]**2)
			else:
				hess[i,j] = (fcn2(h[i],h[j],i,j) - fcn2(h[i],-h[j],i,j) - fcn2(-h[i],h[j],i,j) + fcn2(-h[i],-h[j],i,j))/(4*h[i]*h[j])
				hess[j,i] = hess[i,j]
				
	return hess




