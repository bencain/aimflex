function aim_image_deviate, scaled_pars, $
                            dataimg=dataimg, bkg=bkg, window=window, $
                            psf=psf, fixed_error=fixed_error,$
                            error_scale=error_scale,$
                            sersic=sersic, pseudogauss=pseudogauss,$
                            moffat=moffat,$
                            pranges=pranges, pmids=pmids, $
                            return_norm=return_norm

; This function returns the deviates, pixel by pixel, from a set of
; model parameters to the data image.  This is to be used with
; MPFIT.PRO to fit for the best parameters to fit the data.  The
; statistic used is chi^2=(dataimg-modelimg)^2/error^2, and
; error=sqrt(modelimg). 

; Pick out the pixels where we actually have data.
good=where(window gt 0d,ngood)
if ngood eq 0 then return,9d9

; Undo the -1 to 1 range scaling
p=scaled_pars*pranges+pmids

; Create the model image
modelimg=aim_mk_image(p,bkg,window,psf,$
                     sersic=keyword_set(sersic),$
                     pseudogauss=keyword_set(pseudogauss),$
                     moffat=keyword_set(moffat))

; Compute the difference and error pixel lists and return the
; chi-squared deviation
if not keyword_set(error_scale) then error_scale=1d

diff=(dataimg-modelimg)[good]
errimg=(sqrt(modelimg))[good]*error_scale
if keyword_set(fixed_error) then begin
   if n_elements(fixed_error) ne n_elements(window) then begin
      errimg[*]=abs(fixed_error[0]) 
   endif else errimg[*]=abs(fixed_error[good])
endif

dev=diff/errimg

; If the RETURN_NORM keyword is set, return the chi-squared value
; instead of the deviate.
if keyword_set(return_norm) then return,total(dev^2,/double)

return, dev

end
