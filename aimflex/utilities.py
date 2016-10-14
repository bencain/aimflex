import numpy as np

def set_gaussian_pars(image,model,weights=1.):
	"""
		A function for creating Gaussian model initial parameters from 
		the data and setting them into the model object.
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
	

def STUV(E,g,F,G):
	"""
		Convert ellipticity and lensing parameters to the combined 
		distortion parametrization.
	"""
	M = 1 + E*np.conj(g)
	S = (E+g)/M
	T = (F - np.conj(E)*G)/np.conj(M)
	U = (F - E*np.conj(F))/M
	V = (G - E*F)/M
	
	return S,T,U,V


def calc_moments(image,xord=0,yord=0,weights=1.):
	"""
		A function for calculating arbitrary xy image moments of 
		a 2D image.
	"""
	dims=image.shape[::-1] # PYTHON SHAPES ARE BACKWARDS
	
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


