function mk_gaussian_psf, FWHM

; This function returns a circular Gaussian PSF image.  The size of
; the image is 5 Gaussian widths in radius.

  FWHM=FWHM[0]

; Convert from FWHM to Gaussian width sigma
  sigma=double(FWHM)/sqrt((8d)*alog(2d))

  radius=(5d)*sigma

  size=(2d)*[radius,radius]+1d

; Make the size an integer
  size=floor(size)

; We want an odd number for size 
  size=size+((size+1) mod 2L)

; Mark the midpoint of the array
  ctr=(size-1)/2

; Set the coordinates
  x=double((lindgen(size) mod size[0])-ctr[0])/sigma
  y=double(lindgen(size)/long(size[0])-ctr[1])/sigma

  psf=exp(-(x^2+y^2)/2d)

  outside=where(x^2+y^2 gt (radius/sigma)^2)

  if outside[0] ne -1 then psf[outside]=0d
  
  psf/=total(psf,/double)

  return, psf

end
