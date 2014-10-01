function rp_image_deviate, scaled_pars, $
                           dataimg=dataimg, bkg=bkg, window=window, $
                           psf=psf, fixed_error=fixed_error,$
                           sersic=sersic, pseudogauss=pseudogauss,$
                           doublegauss=doublegauss, moffat=moffat,$
                           pranges=pranges, pmids=pmids, $
                           return_norm=return_norm

; This function returns the deviates, pixel by pixel, from a set of
; model parameters to the data image.  This is to be used with
; MPFIT.PRO to fit for the best parameters to fit the data.  The
; statistic used is chi^2=(dataimg-modelimg)^2/error^2, and
; error=sqrt(modelimg). 

; Pick out the pixels where we actually have data.
good=where(window gt 0d)

; Undo the -1 to 1 range scaling
p=scaled_pars*pranges+pmids

; Create the model image
modelimg=rp_mk_image(p,bkg,window,psf,$
                     sersic=keyword_set(sersic),$
                     pseudogauss=keyword_set(pseudogauss),$
                     doublegauss=keyword_set(doublegauss),$
                     moffat=keyword_set(moffat))

; Compute the difference and error pixel lists and return the
; chi-squared deviation
diff=(dataimg-modelimg)[good]
if keyword_set(fixed_error) then error=fixed_error else error=(sqrt(modelimg))[good]

dev=diff/error

; If the RETURN_NORM keyword is set, return the chi-squared value
; instead of the deviate.
if keyword_set(return_norm) then return,total(dev^2,/double)

return, dev

end
