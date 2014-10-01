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
     print,x[eq1],xiofx[eq1]
  endif

  return,xiofx
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


function aim_nfw_img, theta_s, c200, redshift, coords, $
                      gaussian_pars, window, $
                      lensmodel=lensmodel, kappa=kappa, $
                      ctr=ctr, polar=polar, src_redshift=src_redshift, $
                      tol=tol,$
                      inc_pos=inc_pos, inc_kappa=inc_kappa, sersic=sersic
  
; This is an update that makes NO APPROXIMATIONS in making a data
; image.  One issue that I need to resolve is whether or not I can
; just go pixel by pixel or if I should include a magnification factor
; matching the whole pixel in the image plane back ot the source
; plane.  I likely will.
;
; The lens model parameters can be caught in keywords.

; Some protection vs bad calls

if not keyword_set(ctr) then ctr=[0d,0d]
if not keyword_set(tol) then tol=1d-6

theta_s=double(theta_s)
c200=double(c200)
redshift=double(redshift)

; Set up the outputs
kappa=0d
nlm=6
if keyword_set(inc_pos) then nlm+=2
if keyword_set(inc_kappa) then nlm++
lensmodel=dblarr(nlm)

; Ok, first we make the lensmodel.
if not keyword_set(polar) then begin
   theta=sqrt(total((coords-ctr)^2,/double))
   phi=atan(coords[1]-ctr[1],coords[0]-ctr[0])
endif else begin
   theta=coords[0]
   phi=coords[1]
endelse

; Set the scales:
e_of_z=sqrt(0.3d*(1d + redshift)^3+0.7d)
H_over_c=70.0*e_of_z/3d5 ; Mpc^-1
; Assume a source redshift of 2 unless otherwise specified
if not keyword_set(src_redshift) then src_redshift=2d
D_LS=cosmo_dist([redshift,src_redshift]) ; Each in Mpc
D_S=cosmo_dist(src_redshift)

delta_c=c200^3/(alog(1d + c200) - c200/(1d + c200))

K_s=1.5d*H_over_c^2*D_LS*D_S*delta_c

x=theta/theta_s
if abs(x - 1d) lt tol/2d then begin
   kappa=K_s/3d
   gamma=K_s*(2d*alog(2d) - 5d/3d)
   fflex=0.4d*K_s/theta_s
   gflex=(2d/15d)*(K_s/theta_s)*(60d*alog(2d) - 47d)
endif else begin
   kappa= K_s*(1d - XI(x))/(x^2 - 1d)
   gamma= K_s*(1d - XI(x) - 2d*(1d - 1d/x^2)*(alog(x/2d) + XI(x)))/(x^2 - 1d)
   fflex=-K_s*(2d*x^2 + 1d - 3d*x^2*XI(x))/(theta_s*x*(x^2 - 1d)^2)
   gflex= K_s*(8d*(x - 1d/x)^2*alog(x/2d) + 3d*(1d - 2d*x^2) + $
               (15d*x^2 - 20d + 8d/x^2)*XI(x))/(theta_s*x*(x^2 - 1d)^2)
endelse

g=(gamma/(1d - kappa))*[cos(2d*phi),sin(2d*phi)]
psi1=(fflex/(4d*(1d - kappa)))*[cos(phi),sin(phi)]
psi3=(gflex/(4d*(1d - kappa)))*[cos(3d*phi),sin(3d*phi)]

lensmodel[nlm-6:nlm-1]=[g,psi1,psi3]

if keyword_set(inc_pos) then begin
   lensmodel[0:1]=coords
endif
if keyword_set(inc_kappa) then begin
   lensmodel[nlm-7]=kappa
endif

; Now we make the image

; Unpack the gaussian parameters
logS0=gaussian_pars[0]
alpha=gaussian_pars[3]
thetac=dcomplex(gaussian_pars[1],gaussian_pars[2])
E=dcomplex(gaussian_pars[4],gaussian_pars[5])
I0=(10d)^logS0/(2d*!dpi*alpha^2)

imdims=double(size(window,/dimensions))
dtheta=dcomplex((dindgen(imdims) mod imdims[0]) - (imdims[0] - 1d)/2d,$
                double(lindgen(imdims)/long(imdims[0])) - (imdims[1] - 1d)/2d)

theta0=theta*dcomplex(cos(phi),sin(phi)) - thetac
theta1=theta0 + dtheta

ph0=theta0/sqrt(theta0*conj(theta0))
ph1=theta1/sqrt(theta1*conj(theta1))

x0=real_part(sqrt(theta0*conj(theta0))/theta_s)
x1=real_part(sqrt(theta1*conj(theta1))/theta_s)

a1=(2d*K_s*theta_s*(alog(x1/2d) + XI(x1))/x1)*ph1
a0=(2d*K_s*theta_s*(alog(x0/2d) + XI(x0))/x0)*ph0

dbeta=(dtheta - (a1 - a0))/(1d - kappa)


; Now we're in the source plane
radius=sqrt(real_part( (1d + E*conj(E))*dbeta*conj(dbeta) - $
                       (dbeta^2*conj(E) + conj(dbeta)^2*E) ))

if keyword_set(sersic) then n=gaussian_pars[6] else n=2d

; There is some floating point underflow issue here (probably
; exp(large negative number).  Got this hack from http://www.dfanning.com/
ce=!except
!except=0
void=check_math()
; Make the image

image=I0*exp(-0.5*(radius/alpha)^n)

; Finish the error crap
floating_point_underflow = 32
status = Check_Math()           ; Get status and reset accumulated math error register.
IF(status AND NOT floating_point_underflow) NE 0 THEN begin
   Message, 'AIM_MK_IMAGE.PRO -- IDL Check_Math() error: ' + StrTrim(status, 2)
endif
!except=ce


return,image*window

end
