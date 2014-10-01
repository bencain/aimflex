; FUNCTION NFW_LENSMODEL, COORDS, Z=Z, C200=C200 THETA_S=THETA_S,$
;                         CTR=CTR, POLAR=POLAR, INCL_POS=INCL_POS,$
;                         FORCE_SCALE=FORCE_SCALE
;  
;   There are some extra helper functions defined here, so I've
;   moved a copy of the calling sequence up here for easy reference.
;   I MUST keep this up to date...
;

FUNCTION K,X
  ; Helper function for convergence.  f(y) in Bacon et al. 2005

  gt1=where(x gt 1d,ngt1,complement=le1,ncomplement=nle1)


  fofx=dblarr(n_elements(x))

  if nle1 ne 0 then fofx[le1] = $
     1d - 2d*atanh(sqrt((1d - x[le1])/$
                        (1d + x[le1])))/sqrt(1d - x[le1]^2)

  if ngt1 ne 0 then fofx[gt1] = $
     1d - 2d*atan(sqrt((x[gt1] - 1d)/$
                       (x[gt1] + 1d)))/sqrt(x[gt1]^2 - 1d) 

  return,fofx
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION S,X
  ; Helper function for shear.  g(x) from Golse & Kneib 2002.

  gt1=where(x gt 1d,ngt1,complement=le1,ncomplement=nle1)


  sofx=dblarr(n_elements(x))

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


  hofx=dblarr(n_elements(x))

  if nle1 ne 0 then hofx[le1] = $
     2d*x[le1]*atanh(sqrt((1d - x[le1])/(1d + x[le1])))/sqrt(1d - x[le1]^2) - 1d/x[le1]

  if ngt1 ne 0 then hofx[gt1] = $
     2d*x[gt1]*atan(sqrt((x[gt1] - 1d)/(x[gt1] + 1d)))/sqrt(x[gt1]^2 - 1d) - 1d/x[gt1]

  return,hofx
END
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


FUNCTION G,X
  ; Helper function for 3-flexion. g(y) from Bacon et al. 2005

  gt1=where(x gt 1d,ngt1,complement=le1,ncomplement=nle1)


  gofx=dblarr(n_elements(x))

  if nle1 ne 0 then gofx[le1] = $
     (8d/x[le1]^3 - 20d/x[le1] + 15d*x[le1])*$
     2d*atanh(sqrt((1d - x[le1])/(1d + x[le1])))/sqrt(1d - x[le1]^2)

  if ngt1 ne 0 then gofx[gt1] = $
     (8d/x[gt1]^3 - 20d/x[gt1] + 15d*x[gt1])*$
     2d*atan(sqrt((x[gt1] - 1d)/(x[gt1] + 1d)))/sqrt(x[gt1]^2 - 1d)

  return,gofx
END
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


FUNCTION NFW_LENSMODEL, COORDS, Z=Z, C200=C200, THETA_S=THETA_S,$
                        CTR=CTR, POLAR=POLAR, INCL_POS=INCL_POS,$
                        FORCE_SCALE=FORCE_SCALE
; Put all the coordinates into radial (r,phi) coords.
if not keyword_set(ctr) then ctr=transpose([0d,0d])
if n_elements(ctr) lt 2 then ctr=transpose([ctr[0],ctr[0]])

; Allow for an overall scaling
if not keyword_set(force_scale) then force_scale=1d

if n_elements(size(coords,/dimensions)) lt 2 then coo=transpose(coords) else coo=coords
dims=size(coo,/dimensions)
if dims[1] lt 2 then for i=0,dims[0]-1 do coo[i,1]=coo[i,0]

if not keyword_set(polar) then begin
   coo[*,0]-=ctr[0]
   coo[*,1]-=ctr[1]
   r=sqrt(total(coo^2,2))
   phi=atan(coo[*,1],coo[*,0])
endif else begin
   r=coo[*,0]
   phi=coo[*,1]
endelse

; See if the NFW model is well defined
n_nfw_par=keyword_set(c200)+keyword_set(theta_s)+keyword_set(z)

; Put everything in terms of theta_e and theta_s
if n_nfw_par ne 3 then begin
   print, 'NFW MODEL UNDERDETERMINED!'
   print, '  --- setting THETA_S = 1, c200=5, z=0.2'
   theta_s=1d
   z=0.2d
   c200=5d
endif

; Put radii in units of the scale radius
x=r/theta_s

; Now we want the scale convergence K_s from the redshift, c200 and
; theta_s.
H_over_c=70.0*sqrt(0.3d*(1d + z)^3+0.7d)/3d5 ; Mpc^-1
; Assume a source redshift of 2
D_LS=cosmo_dist([z,2d]) ; Each in Mpc
D_S=cosmo_dist(2d)
D_L=cosmo_dist(z)

K_s=1.5d*200d*(D_LS*D_L^2/D_S)*H_over_c^2*(c200*(1d + c200)^2)*theta_s


; Calculate the convergence
kappa=2d*K_s*K(x)/(x^2 - 1d)

; Calculate the shear
shear=dblarr(dims[0],2)
shear[*,0]=2d*K_s*(2d*S(x)/x^2 - K(x)/(x^2 - 1d))*cos(2d*phi)
shear[*,1]=2d*K_s*(2d*S(x)/x^2 - K(x)/(x^2 - 1d))*sin(2d*phi)


; Calculate 1-flexion
flex1=dblarr(dims[0],2)
flex1[*,0]=((2d*K_s/theta_s)*(2d*x*K(x) - F(x))/(x^2 - 1d))*cos(phi)
flex1[*,1]=((2d*K_s/theta_s)*(2d*x*K(x) - F(x))/(x^2 - 1d))*sin(phi)

; Calculate 3-flexion
flex3=dblarr(dims[0],2)
flex3[*,0]=(2d*K_s/theta_s)*$
           (8d*alog(x/2d)/x^3 + (3d*(1d - 2d*x^2)/x + G(x))/(x^2 - 1d)^2)*cos(3d*phi)
flex3[*,1]=(2d*K_s/theta_s)*$
           (8d*alog(x/2d)/x^3 + (3d*(1d - 2d*x^2)/x + G(x))/(x^2 - 1d)^2)*sin(3d*phi)

; Convert to reduced versions.
lmodel=dblarr(dims[0],6)
lmodel[*,0]=shear[*,0]/(1d - kappa)
lmodel[*,1]=shear[*,1]/(1d - kappa)
lmodel[*,2]=0.25d*flex1[*,0]/(1d - kappa)
lmodel[*,3]=0.25d*flex1[*,1]/(1d - kappa)
lmodel[*,4]=0.25d*flex3[*,0]/(1d - kappa)
lmodel[*,5]=0.25d*flex3[*,1]/(1d - kappa)

; Do the forced scaling
lmodel*=force_scale

if keyword_set(incl_pos) then begin
   lm=dblarr(dims[0],8)
   lm[*,0:1]=coords
   lm[*,2:7]=lmodel
   lmodel=lm
endif

return,lmodel

end

