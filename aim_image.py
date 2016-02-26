import numpy as np

class AIMImage:
	"""
		This class is for AIM fitting images.
	"""
	
	def __init__(self,images):
		"""
			Return an AIMImage object based off the tuple of images: data, window, and 
			psf.  Initial guess parameters are taken from the data image
		"""
		
		# Assumes (data, window, psf) format
		# - > Saves to images attribute
		# The data and window images are assumed to be the same size/shape
		if isinstance(images,(list,tuple)):
			if len(images) > 2 :
				self.images=tuple(images)[0:3] 
			elif len(images)>1 :
				self.images=tuple(images)+(None,)
			elif len(images)>0 :
				self.images=tuple(images)+(None,None,)
			else:
				print "Bad input images for AIMImage"
				self.images=(None,None,None)
			
		else:
			self.images=(images,None,None,)
					
		self.pars_init=self.get_init_pars() 
		self.pars_model=self.pars_init

	
	def get_init_pars(self):
	
		# XY grids
		nx,ny=self.images[0].shape
		linx=np.linspace(-(nx-1)/2.,(nx-1)/2.,nx)
		liny=np.linspace(-(ny-1)/2.,(ny-1)/2.,ny)
		x,y =np.meshgrid(linx,liny,indexing='ij')
				
		# Apply the window
		if self.images[1] != None:
			Iw=self.images[0]*self.images[1]
		else: # or not
			Iw=self.images[0]

		#
		flux=np.sum(Iw)
		c1=np.sum(Iw*x)/flux
		c2=np.sum(Iw*y)/flux
		q11=np.sum(Iw*x**2)/flux - c1**2
		q22=np.sum(Iw*y**2)/flux - c2**2
		q12=np.sum(Iw*x*y)/flux - c1*c2
		
		# Size and ellipticity
		alpha = np.sqrt(np.sqrt(q11*q22-q12**2))
		e1=(q11-q22)/(q11+q22+2*np.sqrt(q11*q22-q12**2))
		e2=2*q12/(q11+q22+2*np.sqrt(q11*q22-q12**2))
		
		print np.log10(flux),c1,c2,alpha,e1,e2
		
		return 0