function aim_sis_img, theta_e, coords, gaussian_pars, window, $
                     lensmodel=lensmodel, kappa=kappa, $
                     ctr=ctr, polar=polar, $
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

theta_e=double(theta_e)

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

kappa=0.5d*theta_e/theta
gamma=(-0.5d*theta_e/theta)*[cos(2d*phi),sin(2d*phi)]
fflex=(-0.5d*theta_e/theta^2)*[cos(phi),sin(phi)]
gflex=(1.5d*theta_e/theta^2)*[cos(3d*phi),sin(3d*phi)]

g=gamma/(1d - kappa)
psi1=fflex/(4d*(1d - kappa))
psi3=gflex/(4d*(1d - kappa))

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

; For an SIS, the deviation is = theta_e*(theta_hat)
dtheta=dcomplex(dx,dy)
theta0=theta*dcomplex(cos(phi),sin(phi))
thetac=dcomplex(xc,yc)

theta1=theta0+dtheta-thetac
norm0=sqrt(theta0*conj(theta0))
norm1=sqrt(theta1*conj(theta1))

;a1=2d*theta_e*theta1/(2d*norm1 - 1d)
;a0=2d*theta_e*theta0/(2d*norm0 - 1d)
;dbeta=dtheta + a0 - a1

; Redid the inclusion of the 1/(1-kappa) so that it matches.  Also
; fixed a stupid error in theta0. 8/4/10
a1=theta_e*(theta1/norm1)
a0=theta_e*(theta0/norm0)
dbeta=(dtheta + (a0 - a1))/(1d - kappa)

; Now we're in the source plane

radius=sqrt(real_part( (1d + E*conj(E))*dbeta*conj(dbeta) - $
                       (dbeta^2*conj(E) + conj(dbeta)^2*E) ))

if keyword_set(sersic) then n=gaussian_pars[6] else n=2d

image=I0*exp(-0.5*(radius/alpha)^n)

return,image*window

end
