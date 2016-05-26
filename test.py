import aim
import numpy as np
from matplotlib import pyplot
from astropy.io import fits
from astropy.modeling import fitting
from astropy.io import ascii

im = '../aimdata/hlsp_frontier_hst_acs-30mas-selfcal_abell2744_f814w_v1.0-epoch2_drz.fits'
wt = '../aimdata/hlsp_frontier_hst_acs-30mas-selfcal_abell2744_f814w_v1.0-epoch2_wht.fits'
cat= '../aimdata/short.cat'

aim.fit_dataset(im, wt, cat, 'dummy.txt', rscale=3.,
				ntag='NUMBER',xtag='X_IMAGE',ytag='Y_IMAGE',atag='A_IMAGE')


# alpha=4.
# Emag = 0.3
# PA = np.pi/2.
# E=Emag*(np.cos(2.*PA) + 1j*np.sin(2.*PA))
# 
# # t_e = 50. #arcsec
# # rad = 55. #arcsec
# phi = 0.  #radian
# # k=(0.5*t_e/rad)
# # g=-(0.5*t_e/rad)*(np.cos(2.*phi) + 1j*np.sin(2.*phi))/(1.-k)
# # F=-0.25*(0.5*t_e/rad**2)*(np.cos(phi) + 1j*np.sin(phi))/(1.-k)
# # G=0.25*(1.5*t_e/rad**2)*(np.cos(3.*phi) + 1j*np.sin(3.*phi))/(1.-k)
# 
# g=-(0.7)*(np.cos(2.*phi) + 1j*np.sin(2.*phi))
# F=-(0.5)*(np.cos(phi) + 1j*np.sin(phi))
# G=(0.5)*(np.cos(3.*phi) + 1j*np.sin(3.*phi))
# 
# F*= 0.01 #convert to per pixel, 10mas pixel scale
# G*= 0.01 #convert to per pixel, 10mas pixel scale
# 
# print np.abs(E), E
# print np.abs(g), g
# print np.abs(F)*alpha, F
# print np.abs(G)*alpha, G
# 
# sextable = ascii.read('data/test.cat')
# i=125
# print sextable['FLUX_BEST'][i]/(2*np.pi*sextable['A_IMAGE'][i]*sextable['B_IMAGE'][i])
# print sextable['X_IMAGE'].astype(int)[i]
# print sextable['Y_IMAGE'].astype(int)[i]
# 
# ngrid=51
# oned=np.linspace(-(ngrid-1)/2,(ngrid-1)/2,ngrid)
# x,y =np.meshgrid(oned,oned)
# 
# flex=aim.AIMGaussian(logI=1.1,alpha=alpha,
# 								E1=E.real,E2=E.imag,
# 								g1=g.real,g2=g.imag,
# 								F1=F.real,F2=F.imag,
# 								G1=G.real,G2=G.imag)
# 
# Im=flex(x,y)
# 
# # pyplot.imshow(Im)
# # pyplot.show()
# 
# data=Im + np.random.randn(ngrid,ngrid)*0.1
# 
# hdu1 = fits.PrimaryHDU(Im)
# hdu1.writeto('img.fits',clobber=True)
# 
# hdu2 = fits.PrimaryHDU(data)
# hdu2.writeto('dat.fits',clobber=True)
# 
# mod=aim.AIMGaussian() #_NEW()
# 
# aim.set_gaussian_pars(Im,mod)
# 
# fitter=fitting.LevMarLSQFitter()
# 
# fit = fitter(mod,x,y,data,maxiter=1000)
# 
# #print fitter.fit_info
# 
# hdu3 = fits.PrimaryHDU(fit(x,y))
# hdu3.writeto('fit.fits',clobber=True)
# 
# hdu4 = fits.PrimaryHDU(fit(x,y)-data)
# hdu4.writeto('res.fits',clobber=True)
# 
# # print "----------"
# # print flex
# # print mod
# # print fit


