import aim_image as ai
import numpy as np
from matplotlib import pyplot
from astropy.modeling import models


ngrid=101
oned=np.linspace(-(ngrid-1)/2,(ngrid-1)/2,ngrid)
x,y =np.meshgrid(oned,oned,indexing='ij')

xbar=5.
ybar=3.
sigmax=5.
sigmay=1.
pa=np.pi/2.

gauss2d = models.Gaussian2D(amplitude=1.,x_mean=xbar,y_mean=ybar,x_stddev=sigmax,y_stddev=sigmay,theta=pa)

# I=np.exp(-((x-xbar)**2/sigmax**2 + (y-ybar)**2/sigmay**2)/2.)/(2*np.pi*sigmax*sigmay)


Im=gauss2d(x,y)

v=ai.AIMImage(Im)
pyplot.imshow(Im)
pyplot.show()