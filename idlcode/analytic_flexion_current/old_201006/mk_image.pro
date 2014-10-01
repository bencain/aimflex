function mk_image, mpars, bkg, window, psf, nopsf=nopsf, noise=noise, $
                   gaussian_noise=gaussian_noise, seed=seed, floor=floor, $
                   flexionscale=flexionscale, sersic=sersic
  
; This function takes in a set of model parameters and returns a model
; image based on the window size.  The model parameters are, in order:
;
; I0     ==> Peak surface brightness of the unlensed Gaussian scaled
;            to the background level.
; xc     ==> X-center of the unlensed image relative to the window center
; yc     ==> Y-Center of the unlensed image relative to the window center
; A      ==> Semimajor-axis of the unlensed elliptical Gaussian
; eplus  ==> Cosine component of the ellipticity
; ecross ==> Sine component of the ellipticity
  
; If the keyword SERSIC is set, then there will be a parameter for the
;        Sersic index
; index ==> Sersic profile index 
;                  I(R)=I_0 * exp( -(R/(sqrt(2)*alpha))^(1/n)) 

; g1    ==> X-component of the reduced lensing shear
; g2    ==> Y-component of the reduced lensing shear
;
; -- Flexion is scaled down by flexionscale.  Default is 10^3
;
; G1_1  ==> X-component of the reduced lensing first-flexion
; G1_2  ==> Y-component of the reduced lensing first-flexion
; G3_1  ==> X-component of the reduced lensing second-flexion
; G3_2  ==> Y-component of the reduced lensing second-flexion
;
; Note that this takes on the conventions of Schneider and Er 2008.
; G1 and G3 are the complex derivatives of g (G1=D*g, G3=Dg).
;
; The model for the unlensed image is 
;
; I(x,y)=I0*exp(-0.5*((x-xc)^2/m1+(y-yc)^2/m2+2*m3*x*y))
;
;   M1, M2 and M3 are defined below.
;
; The "global" parameters are
; PSF         ==> The input point-spread function from seeing.
; window.     ==> The window function used to put an aperture on the data.
; bkg         ==> The background level, in counts per pixel.

npar=n_elements(mpars)

if not keyword_set(flexionscale) then flexionscale=1d3

; Unpack the model parameters
S0=mpars[0]
ctr=dcomplex(mpars[1],mpars[2])

; Convert to the M_i formulation and the psi fields
p=convert_epars(mpars,/pol_to_mi)
p=convert_flexpars(p,/gn_to_psin)
p[npar-4:npar-1]/=flexionscale

; Reduced shear and flexion
shear=dcomplex(p[npar-6],p[npar-5])
psi1=dcomplex(p[npar-4],p[npar-3])
psi3=dcomplex(p[npar-2],p[npar-1])

; Make the coordinates for the model image in the image plane
dims=size(window,/dimensions)

coords=$
   dcomplex((dindgen(dims) mod dims[0]) - (dims[0]-1d)/(2d),$
            double(lindgen(dims)/long(dims[0])) - (dims[1]-1d)/2d)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Only this one or the one below should be uncommented.
;;; This determines whether the center offset is in the 
;;; image plane or in the source plane.  This one is image,
;;; below is source.
coords-=ctr
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Transform to the source frame coordinates
beta=coords-shear*conj(coords)-conj(psi1)*coords^2-$
     (2d)*psi1*coords*conj(coords)-psi3*conj(coords)^2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; See above
;beta-=ctr   ;; UPDATE: We want image plane center offsets!!
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



; Create the image
x=real_part(beta)
y=imaginary(beta)

radius=sqrt(0.5d*(x^2/p[3]+y^2/p[4]+(2d)*p[5]*x*y))
if keyword_set(sersic) then index=mpars[6] else index=0.5d

image=S0*bkg*exp(-radius^(1d/index))

; Add in the background
image+=abs(bkg)

; Convolve with the seeing PSF
if not keyword_set(nopsf) then image=convolve(image,psf) else $
   psf=mk_gaussian_psf(0.25d)

; Add Poisson noise if the keyword is set
if keyword_set(noise) then begin
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
