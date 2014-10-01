function f_of_x,x

  eq_1=where(x eq 1d,n_eq_1)
  gt_1=where(x gt 1d,n_gt_1)
  lt_1=where(x lt 1d,n_lt_1)

  out=x*0d


  if n_lt_1 gt 0 then out[lt_1]=(1 - acosh(1d/x[lt_1])/sqrt(1d - x[lt_1]^2))/(x[lt_1]^2 - 1d)
  if n_eq_1 gt 0 then out[eq_1]=1d/3d
  if n_gt_1 gt 0 then out[gt_1]=(1 - acos(1d/x[gt_1])/sqrt(x[gt_1]^2 - 1d))/(x[gt_1]^2 - 1d)

  return,out
end

function g_of_x,x

  eq_1=where(x eq 1d,n_eq_1)
  gt_1=where(x gt 1d,n_gt_1)
  lt_1=where(x lt 1d,n_lt_1)

  out=x*0d


  if n_lt_1 gt 0 then out[lt_1]=alog(x[lt_1]/2d) + acosh(1d/x[lt_1])/sqrt(1d - x[lt_1]^2)
  if n_eq_1 gt 0 then out[eq_1]=1d + alog(0.5d)
  if n_gt_1 gt 0 then out[gt_1]=alog(x[gt_1]/2d) + acos(1d/x[gt_1])/sqrt(x[gt_1]^2 - 1d)

  return,out
end

function h_of_x,x

  eq_1=where(x eq 1d,n_eq_1)
  gt_1=where(x gt 1d,n_gt_1)
  lt_1=where(x lt 1d,n_lt_1)

  out=x*0d


  if n_lt_1 gt 0 then out[lt_1]=(alog(x[lt_1]/2d))^2 - (acosh(1d/x[lt_1]))^2
  if n_eq_1 gt 0 then out[eq_1]=(alog(0.5d))^2
  if n_gt_1 gt 0 then out[gt_1]=(alog(x[gt_1]/2d))^2 + (acos(1d/x[gt_1]))^2

  return,out
end

function j_of_x,x

; I've defined j(x) = (d/dx)[ f(x) ]

  eq_1=where(x eq 1d,n_eq_1)
  gt_1=where(x gt 1d,n_gt_1)
  lt_1=where(x lt 1d,n_lt_1)

  out=x*0d


  if n_lt_1 gt 0 then out[lt_1]=(3d*x[lt_1]*acosh(1d/x[lt_1])/sqrt(1d - x[lt_1]^2) - 1d/x[lt_1] - 2d*x[lt_1])/(x[lt_1]^2 - 1d)^2
  if n_eq_1 gt 0 then out[eq_1]=-0.4d
  if n_gt_1 gt 0 then out[gt_1]=(3d*x[gt_1]*acos(1d/x[gt_1])/sqrt(x[gt_1]^2 - 1d) - 1d/x[gt_1] - 2d*x[gt_1])/(x[gt_1]^2 - 1d)^2

  return,out


end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function nfw_lens, theta, pars, Z_weight=Z_weight, alpha_only=alpha_only, pars_delta_c=pars_delta_c

; This function returns the lensing fields for a spherical NFW lens.
; The output is: 
;
;    psi, alpha1, alpha2, kappa, gamma1, gamma2, 
;    fflex1, fflex2, gflex1, gflex2
;
; The input is the position (theta) and the parameters which describe
; the lensing system (pars).  If they are 1-D inputs, then theta[0:1]
; give the position of the object and pars[0:1] give the Einstein
; radius for a source at infinity and the softening radius.  Both
; radii are assumed to be in the same units as theta.  If theta is
; multidimensional, then it is taken to be nobj x 2, where nobj is the
; number of lensed images.

; pars gives the details of the lens, in order: scale radius,
; concentration at r200, halo redshift

; The Z_weight keyword gives the scaling from a source at infinity to
; a finite source distance.

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
  
  center=p[0:1]
  rs=p[2]
  c=p[3]
  z=p[4]

  t1=t1-center[0]
  t2=t2-center[1]
  
  r=sqrt(t1^2 + t2^2)
  phi=atan(t2,t1)
  x=r/rs
  
  if keyword_set(pars_delta_c) then delta_c=c else delta_c=(200d/3d)*c^3/(alog(c+1d) - c/(c+1d))

  ; Fixed flat concordance cosmology
  omega_matter=0.272d
  omega_lambda=0.728d

  e_of_z=sqrt(omega_lambda + omega_matter*(1d + z)^3)

  DLS_over_DS=1d - (1d + z)*redshift_to_angdd(z)/((1d4 + 1d)*redshift_to_angdd(1d4))

  ks=1.5d*delta_c*redshift_to_angdd(z)*e_of_z^2*rs * DLS_over_DS

;  print,ks
;  print,delta_c
;  print,redshift_to_angdd(z)
;  print,e_of_z^2
;  print,rs
;  print,DLS_over_DS
;  read,blar



  if keyword_set(alpha_only) then begin
       out=dblarr(nobj,2)

       alpha = 4d*ks*rs*g_of_x(x)/x*dcomplex(cos(phi),sin(phi))

       out[*,0]=real_part(alpha)*Z
       out[*,1]=imaginary(alpha)*Z
    endif else begin

       out=dblarr(nobj,10)
       
       psi=2d*ks*rs^2*h_of_x(x)

       phase=dcomplex(cos(phi),sin(phi))

       alpha = 4d*ks*rs*g_of_x(x)/x*phase

       kappa=2d*ks*f_of_x(x)

       gamma=2d*ks*(2d*g_of_x(x)/x^2 - f_of_x(x))*phase^2

       fflex=(2d*ks/rs)*j_of_x(x)*phase

       gflex=(2d*ks/rs)*(4d*f_of_x(x)/x - j_of_x(x) - 8d*g_of_x(x)/x^3)*phase^3

       out[*,0]=psi*Z
       out[*,1]=real_part(alpha)*Z
       out[*,2]=imaginary(alpha)*Z
       out[*,3]=kappa*Z
       out[*,4]=real_part(gamma)*Z
       out[*,5]=imaginary(gamma)*Z
       out[*,6]=real_part(fflex)*Z
       out[*,7]=imaginary(fflex)*Z
       out[*,8]=real_part(gflex)*Z
       out[*,9]=imaginary(gflex)*Z

    endelse


  if dim_t eq 1 then out=transpose(out)
  
  return,out
  
end
