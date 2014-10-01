function fit_psfmodel, starcat, ngrid, imgdims, npar=npar,$
                       circ=circ, ellip=ellip

; This function fits a psf model from the star catalog and returns a
; set of parameters for making a psf to be convolved with the model
; image.  NGRID is the number of "pixels" for averaging between the
; stellar images.  IMGDIMS is the size of the full image.  The CIRC
; keyword returns the size of a circular Gaussian in each pixel.  The
; ELLIP keyword returns the size, axis ratio and position angle for an
; elliptical Gaussian in each pixel.  If none of these keywords are
; set, then a delta function psf model is returned.  NPAR is a named
; variable which catches the number of psf model parameters.

; Initially we'll just use shape measurements from SExtractor,
; but eventually we'll want to fit to the actual data images.

; MODEL is an array of dimensions NGRID^2 x NPAR+2, where
; MODEL[*,0] are the grid node X values, MODEL[*,1] are the grid
; node Y values and MODEL[*,2:NPAR+1] are the model parameters at
; each grid node.

if keyword_set(circ) then npar=1 else if keyword_set(ellip) then $
   npar=3 else npar=1

model=dblarr(ngrid^2,npar+2)

gridsize=double(imgdims)/double(ngrid)

model[*,0]=((dindgen(ngrid^2) mod ngrid)+0.5d)*gridsize[0]
model[*,1]=double((lindgen(ngrid^2)/long(ngrid))+0.5d)*gridsize[1]

if starcat[0] eq -1 then begin
   model[*,2]=0.25d
   return,model
endif


if keyword_set(circ) then begin

; This is for a circular Gaussian PSF, uniform across the entire
; field pixel. 
   defaultsize=median(starcat[*,4])

   for i=0,ngrid^2-1 do begin
      local=where((abs(model[i,0]-starcat[*,1]) lt gridsize[0]/2d) and $
                  (abs(model[i,1]-starcat[*,2]) lt gridsize[1]/2d),$
                  nlocal)

      if (nlocal lt 5) then model[i,2]=defaultsize else begin
         model[i,2]=median(starcat[local,4])

      endelse
   endfor

endif else if keyword_set(ellip) then begin
   
; This is for an elliptical Gaussian PSF.  We average the size and the
; ellipticity over each field pixel.

   for i=0,ngrid^2-1 do begin
      local=where((abs(model[i,0]-starcat[*,1]) lt gridsize[0]/2d) and $
                  (abs(model[i,1]-starcat[*,2]) lt gridsize[1]/2d),$
                  nlocal)

      if (nlocal lt 1) then model[i,2:4]=[0.25d,1d,0d] else begin
         fwhm=median(starcat[local,4]*sqrt(starcat[local,5]))*$
              sqrt((8d)*alog(2d))
         e1=(((1d)-starcat[local,5]^2)/((1d)-starcat[local,5]^2))*$
            cos((2d)*starcat[local,6])
         e2=(((1d)-starcat[local,5]^2)/((1d)-starcat[local,5]^2))*$
            sin((2d)*starcat[local,6])

         chi1=mean(e1)
         chi2=mean(e2)

         eps=sqrt(((1d)-sqrt(chi1^2+chi2^2))/((1d)+sqrt(chi1^2+chi2^2)))
         psi=(0.5d)*atan(chi2,chi1)
         model[i,2:4]=[fwhm,eps,psi]
         endelse

   endfor

endif else begin

; Default to a delta function.

   model[*,2]=0.25d
endelse

return,model

end



