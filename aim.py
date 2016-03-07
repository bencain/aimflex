
import numpy as np
import astropy.modeling as am


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

	logI  = am.Parameter(default=1.)
	alpha = am.Parameter(default=1.)
	c1 = 	am.Parameter(default=0.)
	c2 = 	am.Parameter(default=0.)
	E1 = 	am.Parameter(default=0.,min=-0.7,max=0.7)
	E2 = 	am.Parameter(default=0.,min=-0.7,max=0.7)
	g1 = 	am.Parameter(default=0.)
	g2 = 	am.Parameter(default=0.)
	F1 = 	am.Parameter(default=0.)
	F2 = 	am.Parameter(default=0.)
	G1 = 	am.Parameter(default=0.)
	G2 = 	am.Parameter(default=0.)
	
			
	@staticmethod
	def evaluate(x,y,logI,alpha,c1,c2,E1,E2,g1,g2,F1,F2,G1,G2):
	
		# Need to get Gaussian2D parameters from my parameters
		amp = np.power(10.,logI)	# Gaussian amplitude
		emag=np.sqrt(E1**2 + E2**2) # Ellipticity magnitude
		q=(1-emag)/(1+emag) 		# Axis ratio b/a
		pa=0.5*np.arctan2(E2,E1) 	# Position angle
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


def initial_parameters(image):
	"""
		A function for creating Gaussian model initial parameters from the data.
	"""
	
	cts = calc_moments(image)
	ctr = [calc_moments(image,xord=1),
		   calc_moments(image,yord=1)]
	Q11 = calc_moments(image,xord=2)/cts
	Q22 = calc_moments(image,yord=2)/cts
	Q12 = calc_moments(image,xord=1,yord=1)/cts
	
	m1=Q11-Q12**2/Q22
	m2=Q22-Q12**2/Q11
	m3=-Q12/(Q11*Q22-Q12**2)
	
	
	

#####################################################################################
#################	UTILITIES	#####################################################
#####################################################################################
	
def calc_moments(image,xord=0,yord=0):
	"""
		A function for calculating arbitrary xy image moments
	"""
	
	if xord < 0:
		xord=0
	if yord < 0:
		yord=0
	dims=image.shape
	
	x,y =np.meshgrid(np.linspace(-(dims[0]-1)/2,(dims[0]-1)/2,dims[0]),
					 np.linspace(-(dims[1]-1)/2,(dims[1]-1)/2,dims[1]))
	
	if (xord+yord) > 1:
		cts = calc_moments(image,xord=0,yord=0)
		ctr = [calc_moments(image,xord=1,yord=0)/cts,
			   calc_moments(image,xord=0,yord=1)/cts]
	else:
		ctr=[0.,0.]
	
	return np.sum(image*np.power((x-ctr[0]),xord)*np.power((y-ctr[1]),yord))
	
	

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
	if not isinstance(epars,collections.Sequence):
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
	
	conds=[
	if (not ae_to_pol) and pol_to_mi=False,
                  mi_to_ae=False,  pol_to_ae=False,
                  mi_to_pol=False, ae_to_mi=False

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

	else if pol_to_mi:
	   alpha=epars[0]
	   eplus=epars[1]
	   ecross=epars[2]

	   emag=sqrt(eplus^2+ecross^2)

		# This is
		# E = ( Q11 - Q22 + 2iQ12 )/( Q11 + Q22 )
		# m1=alpha^2*sqrt(1d -emag^2)/(1d -eplus)
		# m2=alpha^2*sqrt(1d -emag^2)/(1d +eplus)
		# m3=-ecross/(alpha^2*sqrt(1d -emag^2))

		# E = ( Q11 - Q22 + 2iQ12 )/( Q11 + Q22 +2sqrt( Q11*Q22 - Q12^2 ) )
		m1=alpha**2*(1. - emag**2)/(1. + emag**2 - 2.*eplus)
		m2=alpha**2*(1. - emag**2)/(1. + emag**2 + 2.*eplus)
		m3=-2.*ecross/(alpha**2*(1. -emag**2))

		outpars[ellipse_pars]=[m1,m2,m3]
endif


if keyword_set(pol_to_ae) then begin
   alpha=pars[ellipse_pars[0]]
   eplus=pars[ellipse_pars[1]]
   ecross=pars[ellipse_pars[2]]
   emag=sqrt(eplus^2+ecross^2)

; This is
;  E = ( Q11 - Q22 + 2iQ12 )/( Q11 + Q22 )
;   eps=sqrt((1d -emag)/(1d +emag))
;   A=alpha/sqrt(eps)
;   xi=0.5d*atan(ecross,eplus)

; This is
;  E = ( Q11 - Q22 + 2iQ12 )/( Q11 + Q22 +2sqrt( Q11*Q22 - Q12^2 ) )
   eps=(1d - emag)/(1d + emag)
   A=alpha/sqrt(eps)
   xi=0.5d*atan(ecross,eplus)

   outpars[ellipse_pars]=[A,eps,xi]

endif

if keyword_set(mi_to_pol) then begin
   m1=pars[ellipse_pars[0]]
   m2=pars[ellipse_pars[1]]
   m3=pars[ellipse_pars[2]]

; This is
;  E = ( Q11 - Q22 + 2iQ12 )/( Q11 + Q22 )
;   alpha=(m1*m2/(1d - m1*m2*m3^2))^(0.25d)
;   eplus=(m1-m2)/(m1+m2)
;   ecross=-2d*m1*m2*m3/(m1+m2)

; This is
;  E = ( Q11 - Q22 + 2iQ12 )/( Q11 + Q22 +2sqrt( Q11*Q22 - Q12^2 ) )
   alpha=(m1*m2/(1d - m1*m2*m3^2))^(0.25d)
   eplus=(m1 - m2)/(m1 + m2 + 2d*sqrt(m1*m2 - (m1*m2*m3)^2))
   ecross=-2d*m1*m2*m3/(m1 + m2 + 2d*sqrt(m1*m2 - (m1*m2*m3)^2))

   outpars[ellipse_pars]=[alpha,eplus,ecross]

endif

if keyword_set(ae_to_mi) then $
   outpars=convert_epars(convert_epars(pars,/ae_to_pol),/pol_to_mi)

if keyword_set(mi_to_ae) then $
   outpars=convert_epars(convert_epars(pars,/mi_to_pol),/pol_to_ae)

return, outpars