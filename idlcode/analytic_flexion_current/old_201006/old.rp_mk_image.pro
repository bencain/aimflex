function rp_mk_image, mpars, bkg, window, psf, $
                      nopsf=nopsf, $
                      noise=noise, gaussian_noise=gaussian_noise, $
                      seed=seed, $
                      floor=floor, $
                      sersic=sersic

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
;                  I(R) = I0 * exp( -(R/(sqrt(2)*alpha))^(1/n)) 
;           Note that if the Sersic profile is being used, then N0 is
;           no longer the total number of counts.  Instead, the total
;           is: 
;                  N_tot = N0 * n * 4^n * (2n - 1)!
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

npar=n_elements(mpars)

; Unpack the model parameters
logN0=mpars[0]
ctr=dcomplex(mpars[1],mpars[2])

if keyword_set(sersic) then n=mpars[6] else n=0.5d

alpha=mpars[3]
eplus=mpars[4]
ecross=mpars[5]
emagsq=eplus^2+ecross^2

I0=((10d)^logN0/(2d*!dpi*alpha^2)) ;/(n*(4d)^n*gamma(2d*n))

; Reduced shear and flexion (psiN)
shear=dcomplex(mpars[npar-6],mpars[npar-5])
psi1=dcomplex(mpars[npar-4],mpars[npar-3])
psi3=dcomplex(mpars[npar-2],mpars[npar-1])

; Make the coordinates for the model image in the image plane
dims=size(window,/dimensions)

coords=$
   dcomplex((dindgen(dims) mod dims[0]) - (dims[0]-1d)/(2d),$
            double(lindgen(dims)/long(dims[0])) - (dims[1]-1d)/2d)

; The center offset is in the image plane.  
coords-=ctr

; Transform to the source frame coordinates
beta=coords-shear*conj(coords)-conj(psi1)*coords^2-$
     (2d)*psi1*coords*conj(coords)-psi3*conj(coords)^2

; Create the image
x=real_part(beta)
y=imaginary(beta)

radius=sqrt( ((1d - eplus)*x^2 + (1d + eplus)*y^2 -2d*ecross*x*y)/sqrt(1d - emagsq^2) )

image=I0*exp(-0.5*(radius/alpha)^(1d/n))

; Add in the background
image+=bkg

; Convolve with the seeing PSF
if not keyword_set(nopsf) then image=convolve(image,psf) else $
   psf=mk_gaussian_psf(0.25d)

; Add Poisson noise if the keyword is set
if (keyword_set(noise) or keyword_set(gaussian_noise)) then begin
   if keyword_set(gaussian_noise) then begin
      gaussian_noise=double(gaussian_noise)
      image*=((1d)+gaussian_noise*randomn(seed,dims))
   endif else begin
      image=poidev(image,seed=seed)
   endelse
endif

if keyword_set(floor) then image=double(floor(image,/L64))

return,image*window


end
