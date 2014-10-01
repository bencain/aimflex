function mk_psfmodel, pos, modelpars, npar

; This function takes in the position of an image in the field and
; creates a psf image from the given psf model.  MODELPARS is assumed
; to be an N^2 x NPAR+2 array, with MODELPARS[*,*,0] the x
; gridpoints, MODELPARS[*,*,1] the y gridpoints and
; MODELPARS[*,*,2:NPAR+1] the model parameters for each gridpoint.

if npar eq 1 then begin

   dist=sqrt((pos[0]-modelpars[*,0])^2+(pos[1]-modelpars[*,1])^2)

   nearest=where(dist eq min(dist))

   psf=mk_gaussian_psf(modelpars[nearest,2])

endif else if npar eq 3 then begin

   dist=sqrt((pos[0]-modelpars[*,0])^2+(pos[1]-modelpars[*,1])^2)

   nearest=where(dist eq min(dist))

   pars=modelpars[nearest,*]

   A=pars[0]/sqrt((8d)*alog(2d))
   B=A*pars[1]
   xi=pars[2]

; Set the size of the psf image
   radius=(5d)*A
   size=(2d)*[radius,radius]+1d
   size=floor(size)
   size=size+((size+1) mod 2L)

   ctr=(size-1)/2

; Create a coordinate image
   x=double((lindgen(size) mod size[0])-ctr[0])
   y=transpose(x)

   psf=exp(-(0.5d)*((x*cos(xi)+y*sin(xi))^2/A^2+$
                    (y*cos(xi)-x*sin(xi))^2/B^2))

   out=where(x^2+y^2 gt radius^2,nout)
   if nout gt 0 then psf[out]=0d

   psf/=total(psf,/double)

endif else psf=mk_gaussian_psf(0.25d)

return,psf

end

