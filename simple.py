# A simple 2D test class from astropy.modeling

import numpy as np
import astropy.modeling as am


class simpleton(am.Fittable2DModel):

	inputs = ('x','y')
	outputs = ('z',)
	
	a=am.Parameter(default=1.)
	b=am.Parameter(default=1.)		
	
	@staticmethod
	def evaluate(x, y, a, b):
		return (a-b)*x + (a+b)*y