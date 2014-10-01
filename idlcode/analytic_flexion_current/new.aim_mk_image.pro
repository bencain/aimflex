; Profile definition functions

function aim_profile_gaussian, mpars, beta

; This function takes in a set of model parameters (MPARS) and a
; complex coordinate array (BETA) and returns the surface brightness
; associated with the appropriate elliptical Gaussian.  I assume that
; BETA has been calculated from a center-subtracted image plane
; position (therefore mpars[1:2] are ignored)

I0=(10d)^mpars[0]/(2d*!dpi*mpars[3]^2)
E=dcomplex(mpars[4],mpars[5])

rsq=real_part((1d + E*conj(E))*beta*conj(beta) - beta^2*conj(E) - conj(beta)^2*E)

; There is some floating point underflow issue here (probably
; exp(large negative number).  Got this hack from http://www.dfanning.com/
ce=!except
!except=0
void=check_math()
; Make the image
img=I0*exp(-0.5d*rsq/mpars[3]^2)
; Finish the error crap
floating_point_underflow = 32
status = Check_Math()           ; Get status and reset accumulated math error register.
IF(status AND NOT floating_point_underflow) NE 0 THEN begin
   Message, 'AIM_MK_IMAGE.PRO -- IDL Check_Math() error: ' + StrTrim(status, 2)
endif
!except=ce

return,img


end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function aim_profile_moffat, mpars, beta

; This function takes in a set of model parameters (MPARS) and a
; complex coordinate array (BETA) and returns the surface brightness
; associated with the appropriate elliptical Moffat profile.  I assume
; that  BETA has been calculated from a center-subtracted image plane
; position (therefore mpars[1:2] are ignored)

I0=(10d)^mpars[0]/(2d*!dpi*mpars[3]^2)
E=dcomplex(mpars[4],mpars[5])

rsq=real_part((1d + E*conj(E))*beta*conj(beta) - beta^2*conj(E) - conj(beta)^2*E)

return,I0*(1d + rsq/mpars[3]^2)^(-mpars[6])

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function aim_profile_sersic, mpars, beta

; This function takes in a set of model parameters (MPARS) and a
; complex coordinate array (BETA) and returns the surface brightness
; associated with the appropriate elliptical Sersic.  I assume that
; BETA has been calculated from a center-subtracted image plane
; position (therefore mpars[1:2] are ignored)

I0=(10d)^mpars[0]/(2d*!dpi*mpars[3]^2)
E=dcomplex(mpars[4],mpars[5])

rsq=real_part((1d + E*conj(E))*beta*conj(beta) - beta^2*conj(E) - conj(beta)^2*E)

; There is some floating point underflow issue here (probably
; exp(large negative number).  Got this hack from http://www.dfanning.com/
ce=!except
!except=0
void=check_math()
; Make the image
img=I0*exp(-0.5d*((rsq/mpars[3]^2)^(mpars[6]/2d)))
; Finish the error crap
floating_point_underflow = 32
status = Check_Math()           ; Get status and reset accumulated math error register.
IF(status AND NOT floating_point_underflow) NE 0 THEN begin
   Message, 'AIM_MK_IMAGE.PRO -- IDL Check_Math() error: ' + StrTrim(status, 2)
endif
!except=ce

return,img

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function aim_profile_pseudogauss, mpars, beta

; This function takes in a set of model parameters (MPARS) and a
; complex coordinate array (BETA) and returns the surface brightness
; associated with the appropriate elliptical pseudo-Gaussian.  I
; assume that BETA has been calculated from a center-subtracted image
; plane position (therefore mpars[1:2] are ignored)

I0=(10d)^mpars[0]/(2d*!dpi*mpars[3]^2)
E=dcomplex(mpars[4],mpars[5])

x=0.5d*real_part((1d + E*conj(E))*beta*conj(beta) - beta^2*conj(E) - conj(beta)^2*E)/mpars[3]^2

return,I0/(1d + x + mpars[6]*x^2/2d + mpars[7]*x^3/6d + x^4/24d)

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


function aim_lensing_transformation, mpars, theta

; This function takes in a set of model parameters, 

npar=n_elements(mpars)
theta_c=dcomplex(mpars[1],mpars[2]

; Reduced shear and flexion (psiN)
shear=dcomplex(mpars[npar-6],mpars[npar-5])
psi1=dcomplex(mpars[npar-4],mpars[npar-3])
psi3=dcomplex(mpars[npar-2],mpars[npar-1])

; The center offset is in the image plane.  
theta-=theta_c

; Transform to the source frame coordinates
beta=theta - shear*conj(theta) $
     - conj(psi1)*theta^2 - (2d)*psi1*theta*conj(theta) $
     - psi3*conj(theta)^2

return,beta

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


function aim_mk_image, mpars, bkg, window, psf, $
                       nopsf=nopsf, $
                       noise=noise, gaussian_noise=gaussian_noise, $
                       seed=seed, floor=floor, img_border=img_border,$
                       sersic=sersic, pseudogauss=pseudogauss,$
                       moffat=moffat

; This function takes in a set of model parameters and returns a model
; image based on the window size.  The model parameters are, in order:
;
; logS0 
;        ==> Total counts in the unlensed Gaussian.  The maximum of
;            the unlensed profile is I0 = S0/(2 * PI * alpha^2)
; Xc     ==> X-center of the unlensed image relative to the window center
; Yc     ==> Y-Center of the unlensed image relative to the window center
; alpha  ==> Geometric mean of the semimajor- and semi-minor axes of
;            the unlensed elliptical Gaussian 
; eplus  ==> Cosine component of the ellipticity
; ecross ==> Sine component of the ellipticity
;
; If the keyword SERSIC is set, then there will be a parameter for the
;        Sersic index
; n      ==> Sersic profile index 
;                  I(R) = I0 * exp( -0.5 * (R/alpha)^n) 
;
; If the keyword PSEUDOGAUSS is set, then there will be two parameters,
; K1 and K2, for an approximation to a Gaussian:
;                  I(R) = I0 / (1 + X + K1 * X^2/2 + K2 * X^3/6 + X^4/24)
; with X = R^2/(2 * alpha^2)
;
; If the keyword MOFFAT is set, then a Moffat profile will be used
;                 I(R) = I0 / ( 1 + (R/alpha)^2 )^b
;
; g1    ==> Real part of the reduced lensing shear
; g2    ==> Imaginary part of the reduced lensing shear
;
; psi11  ==> Real part of the reduced lensing first-flexion
; psi12  ==> Imaginary part of the reduced lensing first-flexion
; psi31  ==> Real part of the reduced lensing second-flexion
; psi32  ==> Imaginary part of the reduced lensing second-flexion
;
; Note that this takes on the conventions of Schneider and Er 2008.
; G1 and G3 are the complex derivatives of g (G1=D*g, G3=Dg) and psi1
; and psi3 are the coefficients in the lensing equation.
;
; The "global" parameters are
; PSF         ==> The input point-spread function from seeing.
; window.     ==> The window function used to put an aperture on the data.
; bkg         ==> The background level, in counts per pixel.


; Need to have a border so that the PSF won't fall off the edge
; and create weird edge effects.
if keyword_set(nopsf) then psf=mk_gaussian_psf(0.25d)
if not keyword_set(img_border) then $
   img_border=max([size(psf,/dimensions),2]) else $
      img_border=max([max(long(img_border)),2])

npar=n_elements(mpars)

; Make the coordinates for the model image in the image plane
dims=size(window,/dimensions) +2*[img_border,img_border]

coords=$
   dcomplex((dindgen(dims) mod dims[0]) - (dims[0]-1d)/(2d),$
            double(lindgen(dims)/long(dims[0])) - (dims[1]-1d)/2d)

beta=aim_lensing_transformation(mpars,coords)

profname='aim_profile_gaussian'
if keyword_set(pseudogauss) then profname='aim_profile_pseudogauss'
if keyword_set(sersic) then profname='aim_profile_sersic'
if keyword_set(moffat) then profname='aim_profile_moffat'

image=call_function(profname,mpars,beta)

; Add in the background
image+=bkg

; Convolve with the seeing PSF
if not keyword_set(nopsf) then image=convolve(image,psf)

; Add Poisson noise if the keyword is set
if (keyword_set(noise) or keyword_set(gaussian_noise)) then begin
   if keyword_set(gaussian_noise) then begin
      gaussian_noise=double(gaussian_noise)
      image+=gaussian_noise*randomn(seed,dims)
   endif else begin
      image=poidev(image,seed=seed)
   endelse
endif

if keyword_set(floor) then image=double(floor(image,/L64))

return,cut_stamp(image,(dims-1)/2,window)


end
