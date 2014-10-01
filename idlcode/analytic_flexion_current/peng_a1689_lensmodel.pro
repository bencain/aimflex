function peng_a1689_lensmodel, catfile, kappa=kappa, includexy=includexy,$
                               polarcoords=polarcoords, $
                               flexionscale=flexionscale
  
; This function takes in a catalog of objects and returns the lensing 
; parameters from the Peng et al. combined convergence model.  See
; peng_sigma_vs_R/sigma_to_kappa.txt for more details on where this
; comes from.

if not keyword_set(flexionscale) then flexionscale=1d3

catalog=read_data(catfile,size=size,/quiet)

center=[1979.281d,2122.3345d] ; in pixels

x=catalog[*,1]-center[0]
y=catalog[*,2]-center[1]

r=sqrt(x^2+y^2) ; in pixels
phi=atan(y,x)

; The power law is
;
; R>100 arcsec fit
;   log(K) = 6.1961973 - 1.6587973 * log(R)
; R<100 arcsec fit
;   log(K) = 1.5696715 - 0.5796644 * log(R)
;
; Where R is in arcsec.  This means that the changeover radius is

r_split=72.766981d ; in arcsec

; Change over from pixels to arcsec with the object position

scale=0.049d ;arcsec per pixel

r*=scale

; Convert the fits into powerlaw models [A,n], where the lensing
; potential is psi=Ar^n
fit_hi=[6.1961973d,-1.6587973d]
fit_lo=[1.5696715d,-0.5796644d]

powlaw_hi=dblarr(2)
powlaw_lo=dblarr(2)

powlaw_hi[1]=fit_hi[1]+2d
powlaw_lo[1]=fit_lo[1]+2d

powlaw_hi[0]=(2d)*exp(fit_hi[0])/(fit_hi[1]+2d)^2
powlaw_lo[0]=(2d)*exp(fit_lo[0])/(fit_lo[1]+2d)^2

; Now make the lens parameters
lo=where(r lt r_split,nlo,complement=hi,ncomplement=nhi)

kappa=dblarr(nlo+nhi)
gamma=dcomplexarr(nlo+nhi)
fflex=dcomplexarr(nlo+nhi)
gflex=dcomplexarr(nlo+nhi)

if nlo gt 0 then begin

   kappa[lo]=(0.5d)*powlaw_lo[1]^2*powlaw_lo[0]*r[lo]^(powlaw_lo[1]-2d)

   gamma[lo]=dcomplex((0.5d)*powlaw_lo[1]*(powlaw_lo[1]-2d)*powlaw_lo[0]*$
                      r[lo]^(powlaw_lo[1]-2d)*cos((2d)*phi[lo]),$
                      (0.5d)*powlaw_lo[1]*(powlaw_lo[1]-2d)*powlaw_lo[0]*$
                      r[lo]^(powlaw_lo[1]-2d)*sin((2d)*phi[lo]))

   fflex[lo]=dcomplex((0.5d)*powlaw_lo[1]^2*(powlaw_lo[1]-2d)*powlaw_lo[0]*$
                      r[lo]^(powlaw_lo[1]-3d)*cos(phi[lo]),$
                      (0.5d)*powlaw_lo[1]^2*(powlaw_lo[1]-2d)*powlaw_lo[0]*$
                      r[lo]^(powlaw_lo[1]-3d)*sin(phi[lo]))

   gflex[lo]=dcomplex((0.5d)*powlaw_lo[1]*$
                      (powlaw_lo[1]-2d)*(powlaw_lo[1]-4d)*powlaw_lo[0]*$
                      r[lo]^(powlaw_lo[1]-3d)*cos((3d)*phi[lo]),$
                      (0.5d)*powlaw_lo[1]*$
                      (powlaw_lo[1]-2d)*(powlaw_lo[1]-4d)*powlaw_lo[0]*$
                      r[lo]^(powlaw_lo[1]-3d)*sin((3d)*phi[lo]))

endif 

if nhi gt 0 then begin
   
   kappa[hi]=(0.5d)*powlaw_hi[1]^2*powlaw_hi[0]*r[hi]^(powlaw_hi[1]-2d)

   gamma[hi]=dcomplex((0.5d)*powlaw_hi[1]*(powlaw_hi[1]-2d)*powlaw_hi[0]*$
                      r[hi]^(powlaw_hi[1]-2d)*cos((2d)*phi[hi]),$
                      (0.5d)*powlaw_hi[1]*(powlaw_hi[1]-2d)*powlaw_hi[0]*$
                      r[hi]^(powlaw_hi[1]-2d)*sin((2d)*phi[hi]))

   fflex[hi]=dcomplex((0.5d)*powlaw_hi[1]^2*(powlaw_hi[1]-2d)*powlaw_hi[0]*$
                      r[hi]^(powlaw_hi[1]-3d)*cos(phi[hi]),$
                      (0.5d)*powlaw_hi[1]^2*(powlaw_hi[1]-2d)*powlaw_hi[0]*$
                      r[hi]^(powlaw_hi[1]-3d)*sin(phi[hi]))

   gflex[hi]=dcomplex((0.5d)*powlaw_hi[1]*$
                      (powlaw_hi[1]-2d)*(powlaw_hi[1]-4d)*powlaw_hi[0]*$
                      r[hi]^(powlaw_hi[1]-3d)*cos((3d)*phi[hi]),$
                      (0.5d)*powlaw_hi[1]*$
                      (powlaw_hi[1]-2d)*(powlaw_hi[1]-4d)*powlaw_hi[0]*$
                      r[hi]^(powlaw_hi[1]-3d)*sin((3d)*phi[hi]))

endif

shear=gamma/((1d)-kappa)
G1=flexionscale*(fflex+shear*conj(fflex))/((1d)-kappa)
G3=flexionscale*(gflex+shear*fflex)/((1d)-kappa)

lpars=dblarr(size[0],7)

lpars[*,0]=catalog[*,0]
lpars[*,1]=real_part(shear)
lpars[*,2]=imaginary(shear)
lpars[*,3]=real_part(G1)
lpars[*,4]=imaginary(G1)
lpars[*,5]=real_part(G3)
lpars[*,6]=imaginary(G3)

coords=catalog[*,1:2]

if keyword_set(polarcoords) then begin
   x=coords[*,0]-center[0]
   y=coords[*,1]-center[1]
   r=sqrt(x^2+y^2)
   t=atan(y,x)
endif

if keyword_set(includexy) then begin
   out=dblarr(size[0],9)
   out[*,0]=lpars[*,0]
   out[*,1:2]=coords
   out[*,3:8]=lpars[*,1:6]
endif else out=lpars


return,out

end

