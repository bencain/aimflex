import numpy as np
import astropy.modeling as am

# Doing a simple test to try to figure this out

class minimodel(am.FittableModel):
	inputs = ('x',)
	outputs = ('y',)
	
	m = am.Parameter(default=1.)
	b = am.Parameter(default=0.)
	
	
	@staticmethod
	def evaluate(x,m,b):
		y = m*x+b
		return y
		

