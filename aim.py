
import numpy as np
import astropy.modeling as am
from astropy.io import ascii, fits


class AIMGaussian(am.Fittable2DModel):
	"""
		This class implements an elliptical Gaussian source plane model lensed 
		by shear and flexion.  Parameters included:
			logI	- peak surface brightness
			alpha	- size = sqrt(ab) = a*sqrt(E)
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
	"""

	logI  = am.Parameter(default=1.,max=4.)
	alpha = am.Parameter(default=1.,min=0.)
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
	
			
	@staticmethod
	def evaluate(x,y,logI,alpha,c1,c2,E1,E2,g1,g2,F1,F2,G1,G2):
	
		# Need to get Gaussian2D parameters from my parameters
		if logI > 3.:
			amp=1000.
		else:
			amp = np.power(10.,logI)	# Gaussian amplitude
		emag=np.sqrt(E1**2 + E2**2) 	# Ellipticity magnitude
		q=(1.-emag)/(1.+emag) 			# Axis ratio b/a
		pa=0.5*np.arctan2(E2,E1) 		# Position angle
		
		# Gracefully manage choosing the wrong axis as a vs b
		if q < 0:
			q=-q
			pa=pa+np.pi/2

		a=alpha/np.sqrt(q)
		b=a*q
		
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
		return gmodel(beta.real,beta.imag)


class AIMGaussian_NEW(am.Fittable2DModel):
	"""
		This class implements an elliptical Gaussian source plane model lensed 
		by shear and flexion.  Parameters included:
			logI	- peak surface brightness
			A		- size 
			c1		- image plane x location for beta=0
			c2		- image plane y location for beta=0
			S1		- + polarized shear+ellipticity
			S2		- x polarized shear+ellipticity
			T1		- First 1-flexion x
			T2		- First 1-flexion y
			U1		- Second 1-Flexion x
			U2		- Second 1-Flexion y
			V1		- 3-Flexion x
			V2		- 3-Flexion y
	"""

	logI  = am.Parameter(default=1.)
	A  =	am.Parameter(default=1.)
	c1 = 	am.Parameter(default=0.)
	c2 = 	am.Parameter(default=0.)
	S1 = 	am.Parameter(default=0.)
	S2 = 	am.Parameter(default=0.)
	T1 = 	am.Parameter(default=0.)
	T2 = 	am.Parameter(default=0.)
	U1 = 	am.Parameter(default=0.)
	U2 = 	am.Parameter(default=0.)
	V1 = 	am.Parameter(default=0.)
	V2 = 	am.Parameter(default=0.)
	
			
	@staticmethod
	def evaluate(x,y,logI,A,c1,c2,S1,S2,T1,T2,U1,U2,V1,V2):
	
		# Need to get Gaussian2D parameters from my parameters
		amp = np.power(10.,logI)	# Gaussian amplitude
				
		gmodel=am.models.Gaussian2D(amplitude=amp,x_mean=0.,y_mean=0.,
								    x_stddev=A,y_stddev=A,theta=0.)
	
		beta0=c1+1j*c2
		S=S1+1j*S2
		T=T1+1j*T2
		U=U1+1j*U2
		V=V1+1j*V2
		
		coo = x+1j*y - beta0
		coo_c = np.conj(coo)
		
		beta = coo - S*coo_c - np.conj(T)*coo**2 \
				- 2*U*coo*coo_c - V*coo_c**2 
				
		# Now for the model:
		return gmodel(beta.real,beta.imag)


def set_gaussian_pars(image,model,weight=1.):
	"""
		A function for creating Gaussian model initial parameters from the data and
		setting them into the model object.
	"""
	
	cts = calc_moments(image,weight=weight)
	ctr = [calc_moments(image,xord=1,weight=weight)/cts,
		   calc_moments(image,yord=1,weight=weight)/cts]
	Q11 = calc_moments(image,xord=2,weight=weight)/cts
	Q22 = calc_moments(image,yord=2,weight=weight)/cts
	Q12 = calc_moments(image,xord=1,yord=1,weight=weight)/cts
	
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



def fit_dataset(image, weight, catalog, outfile, rscale=2.,
				ntag='NUMBER',xtag='X_IMAGE',ytag='Y_IMAGE',atag='A_IMAGE'):
	"""
		Using a data image, a weight image, and a Source Extractor-like catalog
		fit an AIM profile to each of the objects in the catalog, and store the output.
		
		Assumes the SEcatalog has at least the NUMBER, X_IMAGE, Y_IMAGE, and A_IMAGE 
		columns.
	"""

	sextable = ascii.read(catalog)
	img_hdu = fits.open(image)
	wht_hdu = fits.open(weight)
	
	catno=sextable[ntag]
	cutradii = (sextable[atag]*rscale).astype(int)
	xpix = sextable[xtag].astype(int)
	ypix = sextable[ytag].astype(int)
	
	
	for i in range(cutradii.size):
		print ("Object %i" % catno[i])
		x,y =np.meshgrid(np.linspace(-cutradii[i],cutradii[i],2*cutradii[i]+1),
						 np.linspace(-cutradii[i],cutradii[i],2*cutradii[i]+1))
		mask=(x**2+y**2<=cutradii[i]**2).astype(int)
		
		stamp = img_hdu[0].data[(ypix[i]-cutradii[i]):(ypix[i]+cutradii[i]+1),\
									(xpix[i]-cutradii[i]):(xpix[i]+cutradii[i]+1)]*mask
									
		weight = wht_hdu[0].data[(ypix[i]-cutradii[i]):(ypix[i]+cutradii[i]+1),\
									 (xpix[i]-cutradii[i]):(xpix[i]+cutradii[i]+1)]*mask
		outname='obj_%04i_output' % catno[i]
		
		fits.PrimaryHDU(stamp).writeto(outname+'_stamp.fits',clobber=True)
		fits.PrimaryHDU(weight).writeto(outname+'_weight.fits',clobber=True)
		
		model=AIMGaussian()
		set_gaussian_pars(stamp,model,weight=weight)
		
		print model
		
		model.alpha.max=0.75*cutradii[i]
		model.c1.max=0.5*cutradii[i]
		model.c2.max=0.5*cutradii[i]
		model.c1.min=-0.5*cutradii[i]
		model.c2.min=-0.5*cutradii[i]
		
		fitter=am.fitting.LevMarLSQFitter()

		fit = fitter(model,x,y,stamp,maxiter=1000,weights=weight)
		
		print fit
		print np.sum(weight*(stamp-fit(x,y))**2)
		print stamp.size
		
		fits.PrimaryHDU(fit(x,y)).writeto(outname+'_fit.fits',clobber=True)
		fits.PrimaryHDU(stamp - fit(x,y)).writeto(outname+'_residual.fits',clobber=True)
		
		
	
#####################################################################################
#################	UTILITIES	#####################################################
#####################################################################################
	
def calc_moments(image,xord=0,yord=0,weight=1.):
	"""
		A function for calculating arbitrary xy image moments
	"""
	dims=image.shape
	
	# Handle weights
	if not type(weight) == np.ndarray:
		w = 1.
	else:
		w = weight/np.sum(weight) # Normalize

	if xord < 0: # No negative orders
		xord=0
	if yord < 0:
		yord=0
	
	x,y =np.meshgrid(np.linspace(-(dims[0]-1)/2,(dims[0]-1)/2,dims[0]),
					 np.linspace(-(dims[1]-1)/2,(dims[1]-1)/2,dims[1]))
	
	if (xord+yord) > 1: # Calc ctr and flux for higher moments
		cts = calc_moments(image,xord=0,yord=0,weight=w)
		ctr = [calc_moments(image,xord=1,yord=0,weight=w)/cts,
			   calc_moments(image,xord=0,yord=1,weight=w)/cts]
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
