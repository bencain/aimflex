pro sersic_test

bkg=85d
window=mk_window(75)
psf=mk_gaussian_psf(3d)

true_gaussian_pars=[5d-1,0d,0d,5d,0d,0d,dblarr(6)]
true_sersic_pars=[5d-1,0d,0d,5d,0d,0d,0.5d,dblarr(6)]

dataimg_gaussian=mk_image(true_gaussian_pars,bkg,window,psf,$
                          /noise,gnoise=4d-3)
dataimg_sersic=mk_image(true_sersic_pars,bkg,window,psf,/sersic,$
                        /noise,gnoise=4d-3)

mwrfits,dataimg_gaussian,'st_gaussdata.fits'
mwrfits,dataimg_sersic,'st_sersicdata.fits'
mwrfits,dataimg_gaussian-dataimg_sersic,'st_datadiff.fits'

gaussian_fit=fit_image(dataimg_gaussian,bkg,window,psf,chisq=chisq_g,ntry=3,$
                       dof=dof_g,fixedpars=[6,7],fixedvals=[0d,0d])
sersic_fit=fit_image(dataimg_sersic,bkg,window,psf,chisq=chisq_s,ntry=3,$
                     dof=dof_s,fixedpars=[7,8],fixedvals=[0d,0d],/sersic)

fitimg_gaussian=mk_image(gaussian_fit,bkg,window,psf)
fitimg_sersic=mk_image(sersic_fit,bkg,window,psf,/sersic)

mwrfits,fitimg_gaussian,'st_gaussfit.fits'
mwrfits,fitimg_sersic,'st_sersicfit.fits'
mwrfits,fitimg_gaussian-fitimg_sersic,'st_fitdiff.fits'

print,'Gaussian Fit Parameters'
print,gaussian_fit
print,'True Gaussian Parameters'
print,true_gaussian_pars
print,'True Gaussian Parameters - Fit Parameters'
print,true_gaussian_pars-gaussian_fit
print,'Gaussian chisq/dof'
print,strcompress(string(chisq_g)+'/'+string(dof_g),/remove_all)

print,'Sersic Fit Parameters'
print,sersic_fit
print,'True Sersic Parameters'
print,true_sersic_pars
print,'True Sersic Parameters - Fit Parameters'
print,true_sersic_pars-sersic_fit
print,'Sersic chisq/dof'
print,strcompress(string(chisq_s)+'/'+string(dof_s),/remove_all)


end

