import numpy as np
import astropy.modeling as am


class AIMGaussian(am.Fittable2DModel):

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
				- 2*F*coo*coo_c - G*coo_c**2 #- beta0 
				
		# Now for the model:
		return gmodel(beta.real,beta.imag)
		
		
		
		
	
	
	
	
	
	
	
	
	
	
# class AIMImage:
# 	"""
# 		This class is for AIM fitting images.
# 	"""
# 	
# 	def __init__(self,images):
# 		"""
# 			Return an AIMImage object based off the tuple of images: data, 
# 			window, and psf.  Initial guess parameters are taken from the 
# 			data image.
# 		"""
# 		
# 		# Assumes (data, window, psf) format
# 		# - > Saves to images attribute
# 		# The data and window images are assumed to be the same size/shape
# 		if isinstance(images,(list,tuple)):
# 			if len(images) > 2 :
# 				self.images=tuple(images)[0:3] 
# 			elif len(images)>1 :
# 				self.images=tuple(images)+(None,)
# 			elif len(images)>0 :
# 				self.images=tuple(images)+(None,None,)
# 			else:
# 				print "Bad input images for AIMImage"
# 				self.images=(None,None,None)
# 			
# 		else:
# 			self.images=(images,None,None,)
# 					
# 		self.pars_init=self.get_init_pars() 
# 		self.pars_model=self.pars_init
# 
# 	
# 	def get_init_pars(self):
# 	
# 		# XY grids
# 		nx,ny=self.images[0].shape
# 		linx=np.linspace(-(nx-1)/2.,(nx-1)/2.,nx)
# 		liny=np.linspace(-(ny-1)/2.,(ny-1)/2.,ny)
# 		x,y =np.meshgrid(linx,liny,indexing='ij')
# 				
# 		# Apply the window
# 		if self.images[1] != None:
# 			Iw=self.images[0]*self.images[1]
# 		else: # or not
# 			Iw=self.images[0]
# 
# 		#
# 		flux=np.sum(Iw)
# 		c1=np.sum(Iw*x)/flux
# 		c2=np.sum(Iw*y)/flux
# 		q11=np.sum(Iw*x**2)/flux - c1**2
# 		q22=np.sum(Iw*y**2)/flux - c2**2
# 		q12=np.sum(Iw*x*y)/flux - c1*c2
# 		
# 		# Size and ellipticity
# 		alpha = np.sqrt(np.sqrt(q11*q22-q12**2))
# 		e1=(q11-q22)/(q11+q22+2*np.sqrt(q11*q22-q12**2))
# 		e2=2*q12/(q11+q22+2*np.sqrt(q11*q22-q12**2))
# 		
# 		print np.log10(flux),c1,c2,alpha,e1,e2
# 		
# 		return 0
# 		
# 		