;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION XI,X, TOL=TOL
; Helper function for NFW lensing from Leonard et al. 2010
  if not keyword_set(tol) then tol=1d-6

  eq1=where(abs(x-1d) lt tol/2d,neq1)

  gt1=where(x ge 1d + tol,ngt1)
  lt1=where(x le 1d - tol,nlt1)

  xiofx=x*0d

  if nlt1 ne 0 then xiofx[lt1]=$
     (2d/sqrt(1d - x[lt1]^2))*atanh(sqrt((1d - x[lt1])/(1d + x[lt1])))

  if ngt1 ne 0 then xiofx[gt1]=$
     (2d/sqrt( x[gt1]^2 - 1d))*atan(sqrt((x[gt1] - 1d)/(x[gt1] - 1d)))

  if neq1 gt 0 then begin
     xiofx[eq1]=0.5d*(XI(1d - 2d*tol) + XI(1d + 2d*tol))
  endif

  return,xiofx
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION NFW_LENSMODEL, COORDS, REDSHIFT, C200, THETA_S,$
                        CTR=CTR, POLAR=POLAR, INC_POS=INC_POS,$
                        INC_KAPPA=INC_KAPPA, FORCE_SCALE=FORCE_SCALE, $
                        SOURCE_REDSHIFT=SOURCE_REDSHIFT, TOL=TOL

; Put all the coordinates into radial (r,phi) coords.
if not keyword_set(ctr) then ctr=transpose([0d,0d])
if n_elements(ctr) lt 2 then ctr=transpose([ctr[0],ctr[0]])

; How close to allow x to 1 before defaulting to the unity values
if not keyword_set(tol) then tol=1d-6

; Allow for an overall scaling
if n_elements(force_scale) eq 0 then force_scale=1d

if n_elements(size(coords,/dimensions)) lt 2 then coo=transpose(coords) else coo=coords
dims=size(coo,/dimensions)

theta_s=double(theta_s)
c200=double(c200)
redshift=double(redshift)

theta=dblarr(dims[0])
phi=dblarr(dims[0])
nlm=6
if keyword_set(inc_pos) then nlm+=2
if keyword_set(inc_kappa) then nlm++
lensmodel=dblarr(dims[0],nlm)

; Go into polar coords if not already there
if not keyword_set(polar) then begin
   for i=0d,dims[0]-1 do begin
      theta[i]=sqrt(total((coo[i,0:1]-ctr)^2,/double))
      phi[i]=atan(coo[i,1]-ctr[1],coo[i,0]-ctr[0])
   endfor
endif else begin
      theta=coo[*,0]
      phi=coo[*,1]
endelse

; Set the scales:
e_of_z=sqrt(0.3d*(1d + redshift)^3+0.7d)
H_over_c=70d*e_of_z/2.99792458d5 ; h_70 Mpc^-1
D_L=cosmo_dist(redshift)
; Assume a source redshift of infinity unless otherwise specified
if not keyword_set(src_redshift) then begin
   src_factor=1d
endif else begin
   D_LS=cosmo_dist([redshift,src_redshift]) ; Each in Mpc
   D_S=cosmo_dist(src_redshift)
   src_factor=D_LS/D_S
endelse

delta_c=c200^3/(alog(1d + c200) - c200/(1d + c200))

K_s=1.5d*H_over_c^2*D_L^2*delta_c*src_factor

x=theta/theta_s

kappa=x*0d
gamma=x*0d
fflex=x*0d
gflex=x*0d


eq1=where(abs(x - 1d) lt tol/2d,neq1,complement=ok,ncomplement=nok)

if neq1 gt 0 then begin
   kappa[eq1]=K_s/3d
   gamma[eq1]=K_s*(2d*alog(2d) - 5d/3d)
   fflex[eq1]=0.4d*K_s/theta_s
   gflex[eq1]=(2d/15d)*(K_s/theta_s)*(60d*alog(2d) - 47d)
endif
if nok gt 0 then begin
   kappa[ok]= K_s*(1d - XI(x[ok]))/(x[ok]^2 - 1d)
   gamma[ok]= K_s*(1d - XI(x[ok]) $
                   - 2d*(1d - 1d/x[ok]^2)*(alog(x[ok]/2d) + XI(x[ok])))/(x[ok]^2 - 1d)
   fflex[ok]=-K_s*(2d*x[ok]^2 + 1d - 3d*x[ok]^2*XI(x[ok]))/(theta_s*x[ok]*(x[ok]^2 - 1d)^2)
   gflex[ok]= K_s*(8d*(x[ok] - 1d/x[ok])^2*alog(x[ok]/2d) + 3d*(1d - 2d*x[ok]^2) + $
                   (15d*x[ok]^2 - 20d + 8d/x[ok]^2)*XI(x[ok]))/(theta_s*x[ok]*(x[ok]^2 - 1d)^2)
endif

g=(gamma/(1d - kappa)) * force_scale
psi1=(fflex/(4d*(1d - kappa))) * force_scale
psi3=(gflex/(4d*(1d - kappa))) * force_scale

lensmodel[*,nlm-6]=g*cos(2d*phi)
lensmodel[*,nlm-5]=g*sin(2d*phi)
lensmodel[*,nlm-4]=fflex*cos(phi)
lensmodel[*,nlm-3]=fflex*sin(phi)
lensmodel[*,nlm-2]=gflex*cos(3d*phi)
lensmodel[*,nlm-1]=gflex*sin(3d*phi)

if keyword_set(inc_pos) then begin
   lensmodel[*,0:1]=coords[*,0:1]
endif
if keyword_set(inc_kappa) then begin
   lensmodel[*,nlm-7]=kappa
endif

return,lensmodel

end
