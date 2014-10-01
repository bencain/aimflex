function aim_mk_image, mpars, bkg, window, psf, $
                       nopsf=nopsf, $
                       noise=noise, gaussian_noise=gaussian_noise, $
                       seed=seed, floor=floor, img_border=img_border,$
                       sersic=sersic, pseudogauss=pseudogauss,$
                       moffat=moffat

; This function takes in a set of model parameters and returns a model
; image based on the window size.  The model parameters are, in order:
;
; log10(N0) 
;        ==> Total counts in the unlensed Gaussian.  The maximum of
;            the unlensed Gaussian I0 = N0/(2 * PI * alpha^2)
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
;           Note that if the Sersic profile is being used, then N0 is
;           no longer the total number of counts.
;
; If the keyword PSEUDOGAUSS is set, then there will be two parameters,
; K1 and K2, for an approximation to a Gaussian:
;                  I(R) = I0 / ( 1 + X + K1 * X^2/2 + K2 * X^3/6)^-1
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
if not keyword_set(img_border) then $
   img_border=max([size(psf,/dimensions),2]) else img_border=max([long(img_border[0]),0])



npar=n_elements(mpars)

; Unpack the model parameters
logN0=mpars[0]
ctr=dcomplex(mpars[1],mpars[2])
alpha=mpars[3]
ellip=dcomplex(mpars[4],mpars[5])

I0=((10d)^logN0/(2d*!dpi*alpha^2))

if keyword_set(sersic) then n=mpars[6] else n=2d
if keyword_set(pseudogauss) then begin
   k1=mpars[6]
   k2=mpars[7]
endif

if keyword_set(moffat) then begin
   b=mpars[6]
   I0*=2d*(b-1d)
endif

; Reduced shear and flexion (psiN)
shear=dcomplex(mpars[npar-6],mpars[npar-5])
psi1=dcomplex(mpars[npar-4],mpars[npar-3])
psi3=dcomplex(mpars[npar-2],mpars[npar-1])

; Make the coordinates for the model image in the image plane
dims=size(window,/dimensions) +2*[img_border,img_border]

coords=$
   dcomplex((dindgen(dims) mod dims[0]) - (dims[0]-1d)/(2d),$
            double(lindgen(dims)/long(dims[0])) - (dims[1]-1d)/2d)

; The center offset is in the image plane.  
coords-=ctr

; Transform to the source frame coordinates
beta=coords-shear*conj(coords)-conj(psi1)*coords^2-$
     (2d)*psi1*coords*conj(coords)-psi3*conj(coords)^2

; Create the image
core=0d  ;1d-3

; This is using 
;   E = ( Q11 - Q22 + 2iQ12 )/( Q11 + Q22 )
;radius=real_part($
;       sqrt( beta*conj(beta) - 0.5d*(beta^2*conj(ellip) + conj(beta)^2*ellip)) $
;                )+core

; This is using 
;   E = ( Q11 - Q22 + 2iQ12 )/( Q11 + Q22 +2sqrt( Q11*Q22 - Q12^2 ) )
radius=real_part($
       sqrt( (1d + ellip*conj(ellip))*beta*conj(beta) - $
             (beta^2*conj(ellip) + conj(beta)^2*ellip)$
           )/sqrt(1d - ellip*conj(ellip))$
                )+core


if keyword_set(pseudogauss) then begin
   u=radius^2/(2d*alpha^2)
   image=I0/(1d + u + $
             k1*u^2/2d + $
             k2*u^3/6d + $
             u^4/24d + $
             0d)
endif else if keyword_set(moffat) then begin
   image=I0/(1d + (radius/alpha)^2)^b
endif else begin
; There is some floating point underflow issue here (probably
; exp(large negative number).  Got this hack from http://www.dfanning.com/
   ce=!except
   !except=0
   void=check_math()
; Make the image
   image=I0*exp(-0.5*(radius/alpha)^(n))
; Finish the error crap
   floating_point_underflow = 32
   status = Check_Math()         ; Get status and reset accumulated math error register.
   IF(status AND NOT floating_point_underflow) NE 0 THEN begin
      Message, 'AIM_MK_IMAGE.PRO -- IDL Check_Math() error: ' + StrTrim(status, 2)
   endif
   !except=ce
endelse

; Convolve with the seeing PSF
if not keyword_set(nopsf) then image=convolve(image,psf) else $
   psf=mk_gaussian_psf(0.25d)

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

; Cut out the border and add in the background
image=cut_stamp(image,(dims-1)/2,window)
image+=bkg*window

return,image


end
