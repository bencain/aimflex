function nie_lens, theta, pars, Z_weight=Z_weight, alpha_only=alpha_only

; This function returns the lensing fields for a non-singular
; isothermal ellipse.  The output is:
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
; convergence offset, the core radius, the ellipticity, and the
; position angle (in radians) and the lens weight for a source at infinity.

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
  
  if np lt 8 then stop,'AAAH!  Too few parameters for nie_lens!'


  center=p[0:1]

  t1 = t1-center[0]
  t2 = t2-center[1]


  sigma=p[2]                    ; Einstein radius => 4pi(sigma_v/c)^2, for s=>0
  k0   =p[3]                    ; Mass sheet offset
  s    =p[4]                    ; Core radius
  eps  =p[5]                    ; Ellipticity => (a-b)/(a+b)
  phi  =p[6]                    ; Position angle
  DLI_DI=p[n_elements(p)-1]     ; Lens scaling for src at infinity

  if eps lt 0 then begin
     eps=-eps
     phi+=!dpi/2d
  endif
  if eps gt 0.95d then eps=0.95d

  if eps lt 1d-4 then return,nis_lens(theta, pars, Z_weight=Z, alpha_only=keyword_set(alpha_only))

; Lens parameters
  q = (1d - eps)/(1d + eps)     ; axis ratio
  if q gt (1d - 1d-5) then q = 1d - 1d-5 else if q lt 1d-5 then q=1d-5
  b = sqrt(q)*sigma^2*0.4806d*DLI_DI


; Rotate to a coordinate system where the axes align with the
; ellipse.  We'll calculate the lensing fields in this
; coordinate system and then un-rotate the lensing fields.
  unrotate=dcomplex(cos(phi),sin(phi))

  x = t1*cos(phi) + t2*sin(phi)
  y = t2*cos(phi) - t1*sin(phi)

  w = sqrt(q^2*(s^2 + x^2) + y^2)
  v = sqrt((w + s)^2 + (1d - q^2)*x^2)


  if keyword_set(alpha_only) then begin

     out=dblarr(nobj,2)

     alpha1 = (b/sqrt(1d - q^2)) * atan( x*sqrt(1d - q^2)/(w + s) )
     alpha2 = (b/sqrt(1d - q^2)) * atanh( y*sqrt(1d - q^2)/(w + q^2*s) )

     alpha1+=k0*x
     alpha2+=k0*y

     out[*,0]=(alpha1*cos(phi) - alpha2*sin(phi))*Z
     out[*,1]=(alpha2*cos(phi) + alpha1*sin(phi))*Z

  endif else begin
     
     out=dblarr(nobj,10)

     alpha1 = (b/sqrt(1d - q^2)) * atan( x*sqrt(1d - q^2)/(w + s) )
     alpha2 = (b/sqrt(1d - q^2)) * atanh( y*sqrt(1d - q^2)/(w + q^2*s) )
     
     psi=x*alpha1 + y*alpha2 - b*s*alog(v)

     kappa=0.5d*b/w

     gamma1=-kappa*(x^2 - y^2 + s^2*(1d - q^2))/v^2
     gamma2=-kappa*(2d*x*y)/v^2

     F1=-0.5d*b*q^2*x/w^3
     F2=-0.5d*b*y/w^3

     G1 = ( 2d*x*y*F2 - (x*x - y*y + (1d - q*q)*s*s)*F1 + $
            2d*y*(1d + s/w)*gamma2 - 2d*x*(1d + q*q*s/w)*gamma1)/(v*v) 
     
     G2 = (-(x*x - y*y + (1d - q*q)*s*s)*F2 - 2d*x*y*F1 - $
           2d*x*(1d + q*q*s/w)*gamma2 - 2.0*y*(1d + s/w)*gamma1)/(v*v) 
  

; Add the constant kappa offset terms
     psi+=0.5d*k0*(x^2+y^2)
     alpha1+=k0*x
     alpha2+=k0*y
     kappa+=k0


     out[*,0]=psi

     out[*,1]=(alpha1*cos(1d*phi) - alpha2*sin(1d*phi))*Z
     out[*,2]=(alpha2*cos(1d*phi) + alpha1*sin(1d*phi))*Z

     out[*,3]=kappa*Z

     out[*,4]=(gamma1*cos(2d*phi) - gamma2*sin(2d*phi))*Z
     out[*,5]=(gamma2*cos(2d*phi) + gamma1*sin(2d*phi))*Z

     out[*,6]=(F1*cos(1d*phi) - F2*sin(1d*phi))*Z
     out[*,7]=(F2*cos(1d*phi) + F1*sin(1d*phi))*Z

     out[*,8]=(G1*cos(3d*phi) - G2*sin(3d*phi))*Z
     out[*,9]=(G2*cos(3d*phi) + G1*sin(3d*phi))*Z

  endelse

  if dim_t eq 1 then out=transpose(out)
  
  return,out

end
