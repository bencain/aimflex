function nis_lens, theta, pars, Z_weight=Z_weight, alpha_only=alpha_only

; This function returns the lensing fields for a non-singular
; isothermal sphere.  The output is:
;
;    psi, alpha1, alpha2, kappa, gamma1, gamma2, 
;    fflex1, fflex2, gflex1, gflex2
;
; The input is the position (theta) and the parameters which describe
; the lensing system (pars).  If theta is a 1-D input, then theta[0:1]
; gives the position of the object. If theta is multidimensional, then
; it is taken to be nobj x 2, where nobj is the number of image
; positions. 

; pars gives the details of the lens, in order: center position, the
; velocity dispersion (in 10^3 km/s) of the NIE, a constant
; convergence offset, the core radius and the lens weight for a source at infinity.

; The Z_weight keyword gives the scaling from a source at infinity to
; a finite source distance for each individual image family.

  dim_t=size(theta,/n_dimensions)

  if dim_t eq 1 then begin
     nobj=1
     t1=theta[0]
     t2=theta[1]
  endif else begin
     nobj=n_elements(theta[*,0])
     t1=theta[0:nobj-1,0]
     t2=theta[0:nobj-1,1]
  endelse

  if not keyword_set(Z_weight) then begin
     Z=1d 
  endif else begin
     if n_elements(Z_weight) lt nobj then begin
        Z=Z_weight[0] 
     endif else begin
        Z=Z_weight[0:nobj-1]
     endelse
  endelse

  p=pars
  np=n_elements(p)
  
  if np lt 6 then stop,'AAAH!  Too few parameters for nie_lens!'


  center=p[0:1]

  x = t1-center[0]
  y = t2-center[1]


  sigma=p[2]                    ; Einstein radius => 4pi(sigma_v/c)^2, for s=>0
  k0   =p[3]                    ; Mass sheet offset
  s    =p[4]                    ; Core radius
  DLI_DI=p[n_elements(p)-1]     ; Lens scaling for src at infinity


; Lens parameters
  b = sigma^2*0.4806d*DLI_DI


  w = sqrt(s^2 + x^2 + y^2)
  v = w + s


  if keyword_set(alpha_only) then begin

     out=dblarr(nobj,2)

     alpha1 = b*x/v
     alpha2 = b*y/v

     alpha1+=k0*x
     alpha2+=k0*y

     out[*,0]=alpha1*Z
     out[*,1]=alpha2*Z

  endif else begin
     
     out=dblarr(nobj,10)

     alpha1 = b*x/v
     alpha2 = b*y/v
     
     psi=x*alpha1 + y*alpha2 - b*s*alog(v)

     kappa=0.5d*b/w

     gamma1=-kappa*(x^2 - y^2)/v^2
     gamma2=-kappa*(2d*x*y)/v^2

     F1=-0.5d*b*x/w^3
     F2=-0.5d*b*y/w^3

     G1=((x^2 - y^2)*(2d*x*kappa - F1*v*w) + (2d*x*y)*(F2*v*w - 2d*y*kappa))/(v^3*w)
     G2=((x^2 - y^2)*(2d*y*kappa - F2*v*w) + (2d*x*y)*(2d*x*kappa - F1*v*w))/(v^3*w)

; Add the constant kappa offset terms
     psi+=0.5d*k0*(x^2+y^2)
     alpha1+=k0*x
     alpha2+=k0*y
     kappa+=k0


     out[*,0]=psi

     out[*,1]=alpha1*Z
     out[*,2]=alpha2*Z

     out[*,3]=kappa*Z

     out[*,4]=gamma1*Z
     out[*,5]=gamma2*Z

     out[*,6]=F1*Z
     out[*,7]=F2*Z

     out[*,8]=G1*Z
     out[*,9]=G2*Z

  endelse

  if dim_t eq 1 then out=transpose(out)
  
  return,out

end
