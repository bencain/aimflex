import aim
import numpy as np
from matplotlib import pyplot
from astropy.io import fits
from astropy.modeling import fitting

ngrid=101
oned=np.linspace(-(ngrid-1)/2,(ngrid-1)/2,ngrid)
x,y =np.meshgrid(oned,oned)


flex=aim.AIMGaussian(logI=1.1,alpha=16.,E1=-0.6)#,F1=-0.05,G1=0.05)

Im=flex(x,y)

# pyplot.imshow(Im)
# pyplot.show()

data=Im + np.random.randn(ngrid,ngrid)*0.1

hdu1 = fits.PrimaryHDU(Im)
hdu1.writeto('img.fits',clobber=True)

hdu2 = fits.PrimaryHDU(data)
hdu2.writeto('dat.fits',clobber=True)

mod=aim.AIMGaussian()

aim.set_gaussian_pars(Im,mod)


fitter=fitting.LevMarLSQFitter()

fit = fitter(mod,x,y,data)



hdu3 = fits.PrimaryHDU(fit(x,y))
hdu3.writeto('fit.fits',clobber=True)

hdu4 = fits.PrimaryHDU(fit(x,y)-data)
hdu4.writeto('res.fits',clobber=True)

print "----------"
print flex
print mod
print fit