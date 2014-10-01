function image_deviate, mpars, dataimg=dataimg, bkg=bkg, window=window, $
                        psf=psf, flexionscale=flexionscale,sersic=sersic

; This function returns the deviates, pixel by pixel, from a set of
; model parameters to the data image.  This is to be used with
; MPFIT.PRO to fit for the best parameters to fit the data.  The
; statistic used is chi^2=(dataimg-modelimg)^2/error^2, and
; error=sqrt(modelimg). 

; Pick out the pixels where we actually have data.
good=where(window gt 0d)

p=mpars

; Create the model image
if not keyword_set(flexionscale) then flexionscale=1d3

modelimg=mk_image(p,bkg,window,psf,flexionscale=flexionscale,$
                  sersic=keyword_set(sersic))
   
; Compute the difference and error pixel lists and return the
; chi-squared deviation
diff=(dataimg-modelimg)[good]
error=(sqrt(modelimg))[good]

dev=diff/error

return, dev

end
