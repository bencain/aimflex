pro mk_cluster_image, img_dims, theta_e, ngal, center=center, tag=tag,$
                      cutflex=cutflex, winscale=winscale, subctr=subctr,$
                      subratio=subratio, nonoise=nonoise, gnoise=gnoise,$
                      galsnr=galsnr, nstar=nstar, starsnr=starsnr, $
                      seeing=seeing, bkg=bkg

; This function creates an image of cluster using an SIS lens model
; centered at the center of the image (unless the CENTER keyword is
; otherwise specified).  Galaxies are arranged in the most nearly
; square grid possible, with the galaxies located within 0.5*theta_e
; of the lens center removed.  See below for the assumptions on the
; intrinsic galaxy shapes and magnitudes.

img_dims=long(img_dims)

img_dims+=(img_dims+1L) mod 2L


if not keyword_set(center) then center=[1d,1d]*(img_dims-1L)/2L
if not keyword_set(subctr) then subctr=[0.4d,0.4d]*(img_dims-1L)/2L
if not keyword_set(galsnr) then galsnr=10d
if not keyword_set(nstar) then nstar=10
if not keyword_set(starsnr) then starsnr=galsnr*5d

if not keyword_set(tag) then tag='' else tag+='_'

if not keyword_set(cutflex) then cutflex=0.25d 
; Maximum shear differential across the galaxy.  Galaxies for which
; a*G_n is greater than this will be cut out.

if not keyword_set(winscale) then winscale=25d
                                ; Diameter of the stamp window for
                                ; each galaxy, in units of A.

if not keyword_set(subratio) then subratio=0d
theta_e_sub=subratio*theta_e

; Find the x-y coordinates of the galaxy grid, as well as the polar
; coordinates of the galaxy locations.
ng_side=ceil(sqrt(ngal))

galx=double(floor(((dindgen(ngal) mod ng_side)+0.5d)*$
                  double(img_dims[0])/double(ng_side)))
galy=double(floor(((lindgen(ngal)/long(ng_side))+0.5d)*$
                  double(img_dims[1])/double(ng_side)))

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

; Assume circular galaxies of a given range of sizes and magnitudes
A=(8d)*randomu(seed,ngal)+10d
eps=dblarr(ngal)+1d
xi=dblarr(ngal)

;win=mk_window(winscale*max(A))
side=ceil(winscale*max(A))
side=(side+1L)-(side mod 2L)
win=dblarr(side,side)+1d

if not keyword_set(bkg) then bkg=1000d

color=(2d)*randomu(seed,ngal,/double)

S0=galsnr*sqrt(bkg)*((1d)+randomu(seed,ngal,/double))/(bkg*2d)

; Insert stars
if not keyword_set(seeing) then seeing=3d ; Seeing FWHM

starx=randomu(seed,nstar)*img_dims[0]
stary=randomu(seed,nstar)*img_dims[1]

color_star=(2d)*(randomu(seed,nstar,/double)-0.5d)

S0_star=starsnr*sqrt(bkg)*((1d)+randomu(seed,nstar,/double))/(bkg*2d)

psf=mk_gaussian_psf(seeing)


; Make the images
fullimg_red=dblarr(img_dims)+bkg
fullimg_blue=dblarr(img_dims)+bkg
pars=dblarr(ngal,12)

;Find all the lensing parameters
for i=0,ngal-1 do begin
   ps=convert_epars([S0[i],0d,0d,A[i],eps[i],xi[i],dblarr(6)],/ae_to_pol)

   lps=sis_lensmodel([galrad[i],galphi[i]],theta_e=theta_e,/polar,kappa=k1)

   sub_lps=sis_lensmodel([subrad[i],subphi[i]],theta_e=theta_e_sub,/polar,$
                         kappa=k2)

   ps[6:11]=lens_superposition(lps,sub_lps,k1,k2)

   pars[i,*]=ps
   
endfor

; Find where the images aren't too lensed
good=where((sqrt(pars[*,8]^2+pars[*,9]^2)*A lt cutflex) and $
           (sqrt(pars[*,10]^2+pars[*,11]^2)*A lt cutflex),ngal)

catalog=dblarr(ngal+nstar,15)

for i=0,ngal-1 do begin
   forloop_status,i,ngal,wnum,label='Making Galaxy Stamps'

   catalog[i,0]=i
   catalog[i,1:2]=[galx[good[i]],galy[good[i]]]
   catalog[i,3:14]=pars[good[i],0:11]

   redstamp=mk_image(pars[good[i],*],bkg,win,psf)-bkg
   bluestamp=mk_image([S0[i]*(10d)^(-color[i]/2.5d),$
                       transpose(pars[good[i],1:11])],bkg,win,psf)-bkg

   fullimg_red=insert_stamp(fullimg_red,redstamp,$
                            [galx[good[i]],galy[good[i]]],/sum)
   fullimg_blue=insert_stamp(fullimg_blue,bluestamp,$
                             [galx[good[i]],galy[good[i]]],/sum)
endfor

for i=ngal,ngal+nstar-1 do begin
   forloop_status,i-ngal,nstar,wnum,label='Making Star Stamps'

   catalog[i,0:2]=[i,starx[i-ngal],stary[i-ngal]]
   catalog[i,3:14]=convert_epars([S0_star[i-ngal],0d,0d,$
                                  seeing/sqrt((8d)*alog(2d)),$
                                  1d,0d,dblarr(6)],/ae_to_pol)
   
   redstamp=mk_gaussian_psf(seeing)
   redstamp*=S0_star[i-ngal]/max(redstamp)
   bluestamp=redstamp*(10d)^(-color_star[i-ngal]/2.5d)
   fullimg_red=insert_stamp(fullimg_red,redstamp,$
                            [starx[i-ngal],stary[i-ngal]],/sum)
   fullimg_blue=insert_stamp(fullimg_blue,bluestamp,$
                             [starx[i-ngal],stary[i-ngal]],/sum)
endfor

forloop_status,0,0,wnum,/delete

if keyword_set(noise) then begin
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


   
