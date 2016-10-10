import aimflex
# import numpy as np
# from matplotlib import pyplot
# from astropy.io import fits
# from astropy.modeling import fitting
# from astropy.io import ascii

im = '../../flexion_development/aimdata/hlsp_frontier_hst_acs-30mas-selfcal_abell2744_f814w_v1.0-epoch2_drz.fits'
wt = '../../flexion_development/aimdata/hlsp_frontier_hst_acs-30mas-selfcal_abell2744_f814w_v1.0-epoch2_wht.fits'
cat= '../../flexion_development/aimdata/test.cat'
# cat= 'aimdata/test.cat'

aimflex.fit_dataset(im, wt, cat, 'test_output.txt', rscale=3.,
					ntag='NUMBER',xtag='X_IMAGE',ytag='Y_IMAGE',atag='A_IMAGE',
					outdir='../../flexion_development/out/',save_fig=True,simplex=False)

# 
# alpha=4.
# Emag = 0.3
# PA = np.pi/2.
# E=Emag*(np.cos(2.*PA) + 1j*np.sin(2.*PA))
# 
# t_e = 50. #arcsec
# rad = 55. #arcsec
# phi = 0.  #radian
# k=(0.5*t_e/rad)
# g=-(0.5*t_e/rad)*(np.cos(2.*phi) + 1j*np.sin(2.*phi))/(1.-k)
# F=-0.25*(0.5*t_e/rad**2)*(np.cos(phi) + 1j*np.sin(phi))/(1.-k)
# G=0.25*(1.5*t_e/rad**2)*(np.cos(3.*phi) + 1j*np.sin(3.*phi))/(1.-k)
# 
# g=-(0.7)*(np.cos(2.*phi) + 1j*np.sin(2.*phi))
# F=-(0.5)*(np.cos(phi) + 1j*np.sin(phi))
# G=(0.5)*(np.cos(3.*phi) + 1j*np.sin(3.*phi))
# 
# F*= 0.01 #convert to per pixel, 10mas pixel scale
# G*= 0.01 #convert to per pixel, 10mas pixel scale
# 
# 
# 
# ngrid=51
# 
# # Make the xy grid
# nx = ny = ngrid
# onex=np.linspace(-(nx-1)/2,(nx-1)/2,nx)
# oney=np.linspace(-(ny-1)/2,(ny-1)/2,ny)
# x,y =np.meshgrid(onex,oney)
# 
# npsf = 11
# sigpsf = 2.0
# onep=np.linspace(-(npsf-1)/2,(npsf-1)/2,npsf)
# px,py = np.meshgrid(onep,onep)
# psf = np.exp(-0.5*(px**2 + py**2)/sigpsf**2)
# 
# 
# flex=aim.AIM(logI=1.1,alpha=alpha, index=0.75,
# 				E1=E.real,E2=E.imag,
# 				g1=g.real,g2=g.imag,
# 				F1=F.real,F2=F.imag,
# 				G1=G.real,G2=G.imag)
# 
# 
# Im = flex(x,y,None)
# 
# 
# data = Im + np.random.randn(ngrid,ngrid)*0.1
# 
# hdu1 = fits.PrimaryHDU(Im)
# hdu1.writeto('img.fits',clobber=True)
# 
# hdu2 = fits.PrimaryHDU(data)
# hdu2.writeto('dat.fits',clobber=True)
# 
# mod=aim.AIM()
# 
# aim.set_gaussian_pars(Im,mod)
# 
# # fitter=fitting.LevMarLSQFitter()
# fitter = aim.AIMSimplexLSQFitter()
# 
# fit = fitter(mod,x,y,psf,data,maxiter=1000)
# 
# 
# print fitter.fit_info
# 
# hdu3 = fits.PrimaryHDU(fit(x,y))
# hdu3.writeto('fit.fits',clobber=True)
# 
# hdu4 = fits.PrimaryHDU(fit(x,y)-data)
# hdu4.writeto('res.fits',clobber=True)
# 
# print "----------"
# print flex
# print mod
# print fit
# 
# 
