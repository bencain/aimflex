pro mk_cluster_image, img_dims, theta_e, ngal, center=center, tag=tag,$
                      cutflex=cutflex, winscale=winscale, subctr=subctr,$
                      subratio=subratio, nonoise=nonoise, gnoise=gnoise,$
                      galsnr=galsnr, nstar=nstar, starsnr=starsnr, $
                      seeing=seeing, bkg=bkg, rbkg=rbkg, bbkg=bbkg,$
                      flexionscale=flexionscale, showstamp=showstamp,$
                      fake_a1689=fake_a1689, nogrid=nogrid, $
                      surf_bri_test=surf_bri_test

; This function creates an image of cluster using an SIS lens model
; centered at the center of the image (unless the CENTER keyword is
; otherwise specified).  Galaxies are arranged in the most nearly
; square grid possible, with the galaxies located within 0.5*theta_e
; of the lens center removed.  See below for the assumptions on the
; intrinsic galaxy shapes and magnitudes.

if n_elements(img_dims) ne 2 then img_dims=[img_dims[0],img_dims[0]]

if keyword_set(fake_a1689) then begin
   img_dims=[4221L,4231L]
   center=[1979.281d,2122.3345d]
endif

img_dims=long(img_dims)
img_dims+=(img_dims+1L) mod 2L


if not keyword_set(center) then center=[1d,1d]*(img_dims-1L)/2L
if not keyword_set(subctr) then begin
   angle=2d*!dpi*randomu(seed,/double)
   subctr=0.75d*double((img_dims-1L)/2L)*[cos(angle),sin(angle)]+center
endif
if not keyword_set(galsnr) then galsnr=7d
if not keyword_set(nstar) then nstar=10
if not keyword_set(starsnr) then starsnr=1.5d*galsnr

if not keyword_set(tag) then tag='' else tag+='_'

if not keyword_set(flexionscale) then flexionscale=1d3
if not keyword_set(cutflex) then cutflex=0.025d*flexionscale
; Maximum shear differential across the galaxy.  Galaxies for which
; alpha*G_n is greater than this will be cut out.  Consider a physical
; flexion limit of 0.4/" (as in Leonard et al.). At HST
; resolution this is a flexion of 0.02/pixel.  Since a shear change of
; 0.4 across an image is very large, we'll drop the limit.


if not keyword_set(winscale) then winscale=25d
                                ; Diameter of the stamp window for
                                ; each galaxy, in units of A.

; Size of the substructure SIS
if not keyword_set(subratio) then subratio=0d
theta_e_sub=subratio*theta_e

; Find the x-y coordinates of the galaxy grid, as well as the polar
; coordinates of the galaxy locations.
ng_side=ceil(sqrt(ngal))

galx=double(floor(((dindgen(ngal) mod ng_side)+0.5d)*$
                  double(img_dims[0])/double(ng_side)))
galy=double(floor(((lindgen(ngal)/long(ng_side))+0.5d)*$
                  double(img_dims[1])/double(ng_side)))

if keyword_set(nogrid) then begin
   galx=randomu(seed,ngal,/double)*img_dims[0]
   galy=randomu(seed,ngal,/double)*img_dims[1]
endif

galrad=sqrt((galx-center[0])^2+(galy-center[1])^2)
galphi=atan(galy-center[1],galx-center[0])

subrad=sqrt((galx-subctr[0])^2+(galy-subctr[1])^2)
subphi=atan(galy-subctr[1],galx-subctr[0])


comment='N, X, Y, I0(red), dX, dY, sqrt(AB), eplus, ecross, '+$
        'g1, g2, G1_1, G1_2, G3_1, G3_2.  SIS theta_e='+$
        strcompress(string(theta_e),/remove_all)+' at ('+$
        strcompress(string(center[0])+','+string(center[1]),/remove_all)+$
        ') SIS_sub theta_e='+strcompress(string(theta_e_sub),/remove_all)+$
        ' at ('+strcompress(string(subctr[0])+','+$
                           string(subctr[1]),/remove_all)+')'

; 1 arcsec at HST resolution = 20 pixels (0.05"/px)
; I'll do something a little funky with the ellipticity
; distribution to try to match the A1689 data a bit later.  For
; now we'll just do a uniform in x and E with random E orientation.

A=10d*randomu(seed,ngal)+5d     ; Uniform between 5 and 15
eps=0.2d +randomu(seed,ngal)*0.8d ; Uniform between 0.2 and 1
xi=2d*!dpi*randomu(seed,ngal)

if keyword_set(surf_bri_test) then begin
   A[*]=3d
   eps[*]=1d
   xi[*]=0d
endif


; Do the bkg stuff
if keyword_set(bkg) then begin
   rbkg=bkg
   bbkg=bkg
endif

if not keyword_set(rbkg) then rbkg=1d3
if not keyword_set(bbkg) then bbkg=750d

; Work out the SNR -> S0 transformation

; pick a fiducial size - the median
alpha_med=median(A*sqrt(eps))

; Note that GNOISE is the fractional standard deviation in the image
; as a whole ( gnoise=stddev(img)/mean(img) ).
if not keyword_set(gnoise) then begin
   S0=(galsnr^2/(rbkg*2d*!dpi*alpha_med^2))*$
      (1d +sqrt(1d +8d*alog(2d)*!dpi*alpha_med^2*rbkg/galsnr^2))
endif else begin
   S0=(2d*alog(2d))*(galsnr*gnoise/(1d - galsnr*gnoise))
endelse

S0*=(1d +2d*randomu(seed,ngal))/2d
   
color=(2d)*randomu(seed,ngal,/double)

if keyword_set(surf_bri_test) then begin
   S0=1d-3 + (1d -1d-3)*dindgen(ngal)/double(ngal-1)
   color[*]=0d
endif

; SNR -> S0 for stars
if not keyword_set(seeing) then seeing=3d ; Seeing FWHM
alpha_star=seeing/sqrt((8d)*alog(2d))

starx=randomu(seed,nstar)*img_dims[0]
stary=randomu(seed,nstar)*img_dims[1]

color_star=(2d)*(randomu(seed,nstar,/double)-0.5d)

if not keyword_set(gnoise) then begin
   S0_star=(starsnr^2/(rbkg*2d*!dpi*alpha_star^2))*$
      (1d +sqrt(1d +8d*alog(2d)*!dpi*alpha_star^2*rbkg/starsnr^2))
endif else begin
  S0_star=(2d*alog(2d))*(starsnr*gnoise/(1d - starsnr*gnoise))
endelse

S0_star*=(1d +2d*randomu(seed,nstar))/2d

psf=mk_gaussian_psf(seeing)


; Make the images
fullimg_red=dblarr(img_dims)+rbkg
fullimg_blue=dblarr(img_dims)+bbkg

pars=dblarr(ngal,12)
pars[*,0]=S0
pars[*,3]=A
pars[*,4]=eps
pars[*,5]=xi
for i=0,ngal-1 do pars[i,*]=convert_epars(pars[i,*],/ae_to_pol)

;Find all the lensing parameters   
if keyword_set(fake_a1689) then begin
   lps=peng_powerlaw_new([[galrad],[galphi]],flexionscale=flexionscale,$
                         kappa=k1,/polar)
   lps=lps[*,2:7]

   ; PENG_POWERLAW_NEW automatically includes the postions

endif else begin
   lps=sis_lensmodel([[galrad],[galphi]],theta_e=theta_e,/polar,$
                     kappa=k1,flexionscale=flexionscale)
endelse

sub_lps=sis_lensmodel([[subrad],[subphi]],theta_e=theta_e_sub,/polar,$
                      kappa=k2,flexionscale=flexionscale)

for i=0,ngal-1 do pars[i,6:11]=lens_superposition(lps[i,*],sub_lps[i,*],$
                                                  k1[i],k2[i])

; Find where the images aren't too lensed
good=where((sqrt(pars[*,8]^2+pars[*,9]^2)*pars[*,3] lt cutflex) and $
           (sqrt(pars[*,10]^2+pars[*,11]^2)*pars[*,3] lt cutflex),ngal)


; Find where the images overlap
xya=dblarr(ngal,3)
xya[*,0]=galx[good]
xya[*,1]=galy[good]
xya[*,2]=pars[good,3]

overlap=cat_overlap(xya,factor=2d,noverlap=noverlap)
if overlap[0] ge 0 then good=set_difference(good,good[overlap])
ngal=n_elements(good)

catalog=dblarr(ngal+nstar,15)

for i=0,ngal-1 do begin
   forloop_status,i,ngal,wnum,label='Making Galaxy Stamps'

; Make the window
   side=ceil(winscale*A[good[i]])
   side=(side+1L)-(side mod 2L)
   win=dblarr(side,side)+1d
   
   disp_lo=floor(side/3d)
   disp_hi=2*disp_lo


   catalog[i,0]=i
   catalog[i,1:2]=[galx[good[i]],galy[good[i]]]
   catalog[i,3:14]=pars[good[i],0:11]

   redstamp=mk_image(pars[good[i],*],rbkg,win,psf,$
                     flexionscale=flexionscale)-rbkg
   bluestamp=mk_image([S0[i]*(10d)^(-color[i]/2.5d),$
                       transpose(pars[good[i],1:11])],bbkg,win,psf,$
                      flexionscale=flexionscale)-bbkg

   if keyword_set(showstamp) then $
      disp_scaled_image,redstamp[disp_lo:disp_hi,disp_lo:disp_hi],wnum_stamp

   fullimg_red=insert_stamp(fullimg_red,redstamp,$
                            [galx[good[i]],galy[good[i]]],/sum)
   fullimg_blue=insert_stamp(fullimg_blue,bluestamp,$
                             [galx[good[i]],galy[good[i]]],/sum)
endfor

if keyword_set(showstamp) then disp_scaled_image,0,wnum_stamp,/delete

for i=ngal,ngal+nstar-1 do begin
   forloop_status,i-ngal,nstar,wnum,label='Making Star Stamps'

   catalog[i,0:2]=[i,starx[i-ngal],stary[i-ngal]]
   catalog[i,3:14]=convert_epars([S0_star[i-ngal],0d,0d,$
                                  seeing/sqrt((8d)*alog(2d)),$
                                  1d,0d,dblarr(6)],/ae_to_pol)
   
   redstamp=mk_gaussian_psf(seeing)
   redstamp*=rbkg*S0_star[i-ngal]/max(redstamp)
   bluestamp=redstamp*(10d)^(-color_star[i-ngal]/2.5d)
   fullimg_red=insert_stamp(fullimg_red,redstamp,$
                            [starx[i-ngal],stary[i-ngal]],/sum)
   fullimg_blue=insert_stamp(fullimg_blue,bluestamp,$
                             [starx[i-ngal],stary[i-ngal]],/sum)
endfor

forloop_status,0,0,wnum,/delete

if not keyword_set(nonoise) then begin
   if keyword_set(gnoise) then begin
      gnoise=double(gnoise)
      fullimg_red*=((1d)+gnoise*randomn(seed,img_dims))
      fullimg_blue*=((1d)+gnoise*randomn(seed,img_dims))
   endif else begin
      fullimg_red=poidev(fullimg_red,seed=seed)
      fullimg_blue=poidev(fullimg_blue,seed=seed)
   endelse
endif


mwrfits,fullimg_red,tag+'clus_img_red.fits'
mwrfits,fullimg_blue,tag+'clus_img_blue.fits'

save_data,catalog,tag+'clus_pars.dat',comment=comment

end


   
