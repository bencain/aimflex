import aim_image as ai
import numpy as np
from matplotlib import pyplot
from astropy.io import fits


ngrid=101
oned=np.linspace(-(ngrid-1)/2,(ngrid-1)/2,ngrid)
x,y =np.meshgrid(oned,oned)


flex=ai.AIMGaussian(logI=0.,alpha=6.,c1=0.,c2=0.,g1=-0.6,F1=-0.05,G1=0.15)

Im=flex(x,y)

pyplot.imshow(Im)
pyplot.show()

hdu = fits.PrimaryHDU(Im)

hdu.writeto('new.fits',clobber=True)