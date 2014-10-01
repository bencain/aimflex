;function aim_nfw_img, theta_s, c200, redshift, coords, gaussian_pars, window, $
;                     lensmodel=lensmodel, kappa=kappa, $
;                     ctr=ctr, polar=polar, src_redshift=src_redshift, $
;                     inc_pos=inc_pos, inc_kappa=inc_kappa, sersic=sersic
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


function aim_nfw_img, theta_s, c200, redshift, coords, gaussian_pars, window, $
                     lensmodel=lensmodel, kappa=kappa, $
                     ctr=ctr, polar=polar, src_redshift=src_redshift, $
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
kappa=2d*K_s*K(x)/(x^2 - 1d)
gamma=2d*K_s*(2d*S(x)/x^2 - K(x))
fflex=(-2d*(K_s/theta_s)/(x^2 - 1d)^2)*(2d*x*K(x) - F(x))
gflex=(2d*(K_s/theta_s))*$
      (8d*alog(x/2d)/x^3 + $
       ((3d/x)*(1d - 2d*x^2) + G(x))/(x^2 - 1d)^2)

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
logN0=gaussian_pars[0]
xc=gaussian_pars[1]
yc=gaussian_pars[2]
alpha=gaussian_pars[3]
E=dcomplex(gaussian_pars[4],gaussian_pars[5])
I0=(10d)^logN0/(2d*!dpi*alpha^2)

imdims=double(size(window,/dimensions))
dx=(dindgen(imdims) mod imdims[0]) - (imdims[0] - 1d)/2d - xc
dy=double(lindgen(imdims)/long(imdims[0])) - (imdims[1] - 1d)/2d - yc

dtheta=dcomplex(dx,dy)
theta0=theta*dcomplex(cos(phi),sin(phi))

thetac=dcomplex(xc,yc)
theta1=theta0+dtheta-thetac

norm0=sqrt(theta0*conj(theta0))
norm1=sqrt(theta1*conj(theta1))

;kappa1=2d*K_s*K(norm1/theta_s)/((norm1/theta_s)^2 - 1d)
;kappa0=2d*K_s*K(norm0/theta_s)/((norm0/theta_s)^2 - 1d)
;a1=4d*K_s*theta_s^2*S(norm1/theta_s)*theta1/(norm1^2*(1d - kappa1))
;a0=4d*K_s*theta_s^2*S(norm0/theta_s)*theta0/(norm0^2*(1d - kappa0))
;dbeta=dtheta + a0 - a1

a1=4d*K_s*theta_s^2*S(norm1/theta_s)*theta1/norm1^2
a0=4d*K_s*theta_s^2*S(norm0/theta_s)*theta0/norm0^2
dbeta=(dtheta + (a0 - a1))/(1d - kappa)


; Now we're in the source plane
radius=sqrt(real_part( (1d + E*conj(E))*dbeta*conj(dbeta) - $
                       (dbeta^2*conj(E) + conj(dbeta)^2*E) ))

if keyword_set(sersic) then n=gaussian_pars[6] else n=2d

image=I0*exp(-0.5*(radius/alpha)^n)

return,image*window

end
