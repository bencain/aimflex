;FUNCTION NFW_LENSMODEL, COORDS, REDSHIFT, C200, THETA_S,$
;                        CTR=CTR, POLAR=POLAR, INC_POS=INC_POS,$
;                        INC_KAPPA=INC_KAPPA, FORCE_SCALE=FORCE_SCALE, $
;                        SOURCE_REDSHIFT=SOURCE_REDSHIFT
;
;   There are some extra helper functions defined here, so I've
;   moved a copy of the calling sequence up here for easy reference.
;   I MUST keep this up to date...
;

FUNCTION K,X
  ; Helper function for convergence.  f(y) in Bacon et al. 2005

  gt1=where(x gt 1d,ngt1,complement=le1,ncomplement=nle1)

  kofx=x*0d

  if nle1 ne 0 then kofx[le1] = $
     1d - 2d*atanh(sqrt((1d - x[le1])/$
                        (1d + x[le1])))/sqrt(1d - x[le1]^2)

  if ngt1 ne 0 then kofx[gt1] = $
     1d - 2d*atan(sqrt((x[gt1] - 1d)/$
                       (x[gt1] + 1d)))/sqrt(x[gt1]^2 - 1d) 

  return,kofx
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION S,X
  ; Helper function for shear.  g(x) from Golse & Kneib 2002.

  gt1=where(x gt 1d,ngt1,complement=le1,ncomplement=nle1)


  sofx=x*0d

  if nle1 ne 0 then sofx[le1] = $
     alog(x[le1]/2d) + acosh(1d/x[le1])/sqrt(1d - x[le1]^2)

  if ngt1 ne 0 then sofx[gt1] = $
     alog(x[gt1]/2d) + acos(1d/x[gt1])/sqrt(x[gt1]^2 - 1d)

  return,sofx
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION F,X
  ; Helper function for 1-flexion.  h(y) in Bacon et al. 2005

  gt1=where(x gt 1d,ngt1,complement=le1,ncomplement=nle1)


  fofx=x*0d

  if nle1 ne 0 then fofx[le1] = $
     2d*x[le1]*atanh(sqrt((1d - x[le1])/(1d + x[le1])))/sqrt(1d - x[le1]^2) - 1d/x[le1]

  if ngt1 ne 0 then fofx[gt1] = $
     2d*x[gt1]*atan(sqrt((x[gt1] - 1d)/(x[gt1] + 1d)))/sqrt(x[gt1]^2 - 1d) - 1d/x[gt1]

  return,fofx
END
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


FUNCTION G,X
  ; Helper function for 3-flexion. g(y) from Bacon et al. 2005

  gt1=where(x gt 1d,ngt1,complement=le1,ncomplement=nle1)


  gofx=x*0d

  if nle1 ne 0 then gofx[le1] = $
     (8d/x[le1]^3 - 20d/x[le1] + 15d*x[le1])*$
     2d*atanh(sqrt((1d - x[le1])/(1d + x[le1])))/sqrt(1d - x[le1]^2)

  if ngt1 ne 0 then gofx[gt1] = $
     (8d/x[gt1]^3 - 20d/x[gt1] + 15d*x[gt1])*$
     2d*atan(sqrt((x[gt1] - 1d)/(x[gt1] + 1d)))/sqrt(x[gt1]^2 - 1d)

  return,gofx
END
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



FUNCTION NFW_LENSMODEL, COORDS, REDSHIFT, C200, THETA_S,$
                        CTR=CTR, POLAR=POLAR, INC_POS=INC_POS,$
                        INC_KAPPA=INC_KAPPA, FORCE_SCALE=FORCE_SCALE, $
                        SOURCE_REDSHIFT=SOURCE_REDSHIFT

; Put all the coordinates into radial (r,phi) coords.
if not keyword_set(ctr) then ctr=transpose([0d,0d])
if n_elements(ctr) lt 2 then ctr=transpose([ctr[0],ctr[0]])

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
   for i=0,dims[0]-1 do begin
      theta[i]=sqrt(total((coo[i,0:1]-ctr)^2,/double))
      phi[i]=atan(coo[i,1]-ctr[1],coo[i,0]-ctr[0])
   endfor
endif else begin
   for i=0,dims[0]-1 do begin
      theta[i]=coo[i,0]
      phi[i]=coo[i,1]
   endfor
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
kappa=2d*K_s*K(x)/(x^2 - 1d)
gamma=2d*K_s*(2d*S(x)/x^2 - K(x))
fflex=(-2d*(K_s/theta_s)/(x^2 - 1d)^2)*(2d*x*K(x) - F(x))
gflex=(2d*(K_s/theta_s))*$
      (8d*alog(x/2d)/x^3 + $
       ((3d/x)*(1d - 2d*x^2) + G(x))/(x^2 - 1d)^2)

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


