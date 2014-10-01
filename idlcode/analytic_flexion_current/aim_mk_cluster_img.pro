PRO AIM_MK_CLUSTER_IMG, IMSIZE, TAG, SEED, $
                        NGAL_BKG=NGAL_BKG, NGAL_MEM=NGAL_MEM,$
                        IMG_NOISE=IMG_NOISE, BKG=BKG,$
                        SERSIC=SERSIC, LENSPLANE_SCALE=LENSPLANE_SCALE,$
                        PSFFWHM=PSFFWHM, SUBSTRUCTURE=SUBSTRUCTURE

if n_elements(imsize) lt 2 then imsize=imsize[0]*[1,1]
imsize=long(imsize)
imsize+=(imsize+1) mod 2

if not keyword_set(lensplane_scale) then lensplane_scale=1d
if not keyword_set(ngal_bkg) then $
   ngal_bkg=ceil(50d*product(imsize*lensplane_scale)*(0.05d/60d)^2)
                                ; Assume 50 bkg gals/sq-arcmin unless
                                ; specified, and HST plate scale
if ngal_bkg lt 1 then $
   ngal_bkg=ceil(50d*product(imsize*lensplane_scale)*(0.05d/60d)^2)
if not keyword_set(ngal_mem) then ngal_mem=0
if ngal_mem lt 0 then ngal_mem=0

if not keyword_set(psffwhm) then psffwhm=1.8d
if psffwhm lt 1.8d then psffwhm=1.8d
if not keyword_set(bkg) then bkg=0d

; Cosmological Information
Ds=cosmo_dist(1d)               ; distances in Mpc
Dl=cosmo_dist(0.2d)
Dls=cosmo_dist([0.2d,1d])

cps_per_kappa=1d                ; Kinda kludgy, but we'll see
bkg_color=2d                    ; All bkg galaxies will be redder
mem_color=0d                    ; All mem galaxies will be bluer

; Data saving places
kappaDM=dblarr(imsize)
kappaMEM=dblarr(imsize)

memcat=dblarr(max([ngal_mem,1]),9)
memcom='Member Galaxies: X1, X2, logS0_blue, logS0_red, alpha, E+, Ex, n, theta_e'
bkgcat=dblarr(ngal_bkg,10)
bkgcom='Background Galaxies: X1, X2, B1, B2, logS0_blue, logS0_red, alpha, E+, Ex, n'


; Image plane coordinates
theta=dcomplex(dindgen(imsize) mod imsize[0],$
               double(lindgen(imsize)/long(imsize[0])))


; Dark Matter Distribution
print,'Creating the DM distribution'
theta_e=20d/0.05d ; 30" Einstein radius in HST pix
theta_c=7d/0.05d ; 10" Core radius in HST pix 
dmctr=dcomplex(0.5d*(imsize[0]-1d),0.5d*(imsize[1]-1d))
kappaDM=theta_e/sqrt(abs(theta-dmctr)^2+theta_c^2)
delta=kappa_to_defl(theta-dmctr,kappaDM)

if keyword_set(substructure) then begin
   print,'Adding a substructure'
   r=0.2d
   theta_e*=r
   theta_c*=r ; (r fraction of mass in substructure)
   dmctr=dcomplex(0.5d*(imsize[0]-1d),0.5d*(imsize[1]-1d))
   dmctr+=dcomplex(0.25d*(imsize[0]-1d),0.25d*(imsize[1]-1d))
   kappaDM=(1d -r)*kappaDM + theta_e/sqrt(abs(theta-dmctr)^2+theta_c^2)
   delta+=kappa_to_defl(theta-dmctr,kappaDM)
endif

Imem_b    =dblarr(imsize)
Imem_r    =dblarr(imsize)

; Do the cluster member stuff
if ngal_mem gt 0 then begin
   print,'Adding the cluster member galaxies'
   theta_gal=dcomplex(0.1d*imsize[0]*randomn(seed,ngal_mem),$
                      0.1d*imsize[1]*randomn(seed,ngal_mem))+dmctr
   bad=where(abs(theta_gal-dmctr) gt 0.5d*min(imsize),nbad)
   while nbad gt 0 do begin
      theta_gal[bad]=dcomplex(0.1d*imsize[0]*randomn(seed,nbad),$
                              0.1d*imsize[1]*randomn(seed,nbad))+dmctr
      bad=where(abs(theta_gal-dmctr) gt 0.5d*min(imsize),nbad)
   endwhile

   E=dcomplex(0.3d*randomn(seed,ngal_mem),0.3d*randomn(seed,ngal_mem))
   bad=where(abs(E) gt 0.6d,nbad)
   while nbad gt 0 do begin
      E[bad]=dcomplex(0.3d*randomn(seed,nbad),0.3d*randomn(seed,nbad))
      bad=where(abs(E) gt 0.6d,nbad)
   endwhile
   
   alpha=5d*(abs(randomn(seed,ngal_mem)*2d + 3d)+4d)
   n=randomu(seed,ngal_mem)*2.5d + 0.2d

   sigv=220d + 20d*randomn(seed,ngal_mem)
   theta_e=4d*!dpi*(sigv/3d5)^2*(Dls/Ds) ; In radians
   theta_e*=(180d/!dpi)*(3600d/0.05d)    ; In HST pix

   logS0_r=alog10(cps_per_kappa*theta_e)
   logS0_b=logS0_r-0.4d*mem_color

   print,'For-loop starting '+systime()
   for i=0,ngal_mem-1 do begin
      forloop_status,i,ngal_mem,wn,label='Member Galaxies'
      dt=theta-theta_gal[i]
      r2=real_part((1d + abs(E[i]^2))*abs(dt)^2 - dt^2*conj(E[i]) - conj(dt)^2*E[i])
      
      core=0.5d                 ; Core for the SIS
      kappa_gal=0.5d*theta_e[i]/sqrt(r2 + core^2)
      kappaMEM+=kappa_gal
      delta+=kappa_to_defl(theta-theta_gal[i],kappa_gal)
      Imem_r+=((10d)^logS0_r[i]/(2d*!dpi*alpha[i]^2))*exp(-0.5d*(r2/alpha[i]^2)^(n[i]/2d))
      
   endfor
   Imem_b=Imem_r*(10d)^(-0.4d*mem_color)

   memcat[*,0]=real_part(theta_gal)
   memcat[*,1]=imaginary(theta_gal)
   memcat[*,2]=logS0_b
   memcat[*,3]=logS0_r
   memcat[*,4]=alpha
   memcat[*,5]=real_part(E)
   memcat[*,6]=imaginary(E)
   memcat[*,7]=n
   memcat[*,8]=theta_e
endif

; Make a convergence map and member galaxy images
kappa=kappaDM+kappaMEM
; Mass
mwrfits,kappaMEM,tag+'_kappaMEM.fits'
mwrfits,kappaDM, tag+'_kappaDM.fits'
mwrfits,kappa,   tag+'_kappa.fits'
delvarx,kappaMEM,kappaDM,kappa ; Clear out memory

; Convert the convergence map to a deflection map
print,'Converting the convergence to deflection'
;delta=kappa_to_defl(theta,kappa)
beta=theta-delta

; Now we need a set of source galaxies
print,'Adding the background source galaxies'
Isrc_ul_r =dblarr(imsize)
Isrc_wl_r =dblarr(imsize)
Isrc_ul_b =dblarr(imsize)
Isrc_wl_b =dblarr(imsize)

; Sources could be off the image plane and lensed in
beta_gal=dcomplex(imsize[0]*(lensplane_scale*randomu(seed,ngal_bkg)-(lensplane_scale-1d)/2d),$
                  imsize[1]*(lensplane_scale*randomu(seed,ngal_bkg)-(lensplane_scale-1d)/2d))

E=dcomplex(0.3d*randomn(seed,ngal_bkg),0.3d*randomn(seed,ngal_bkg))
bad=where(abs(E) gt 0.6d,nbad)
while nbad gt 0 do begin
   E[bad]=dcomplex(0.3d*randomn(seed,nbad),0.3d*randomn(seed,nbad))
   bad=where(abs(E) gt 0.6d,nbad)
endwhile

alpha=2d*(abs(randomn(seed,ngal_bkg)*2d + 3d)+4d)
alpha*=(Dl/Ds)
n=randomn(seed,ngal_bkg) + 2d    ; Normal dist about 2 with sigma=1
bad=where((n gt 3d) or (n lt 0.2d),nbad)
while nbad gt 0 do begin
   n[bad]=randomn(seed,nbad) + 2d    ; Normal dist about 2 with sigma=1
   bad=where((n gt 3d) or (n lt 0.2d),nbad)
endwhile

sigv=220d + 20d*randomn(seed,ngal_bkg)
theta_e=4d*!dpi*(sigv/3d5)^2*(Dls/Ds) ; In radians
theta_e*=(180d/!dpi)*(3600d/0.05d)    ; In HST pix

logS0_r=alog10(cps_per_kappa*theta_e*(Dl/Ds)^2)
logS0_b=logS0_r-0.4d*bkg_color
theta_gal=beta_gal*0

print,'For-loop starting '+systime()
for i=0,ngal_bkg-1 do begin
   forloop_status,i,ngal_bkg,wn,label='Background Galaxies'

   dt_ul=theta-beta_gal[i]      ; Position without lensing
   dt_wl= beta-beta_gal[i]      ; Position WITH lensing

   okscl=15d

   mag_dt_wl=abs(dt_wl)
   ok_wl=where(mag_dt_wl lt okscl^(2d/n[i])*alpha[i],nok_wl)

   mag_dt_ul=abs(dt_ul)
   ok_ul=where(mag_dt_ul lt okscl^(2d/n[i])*alpha[i],nok_ul)

   tmp=min(mag_dt_wl,minloc)
   theta_gal=theta[minloc]
  
   if nok_wl gt 0 then begin
      r2_wl=(1d + abs(E[i]^2))*abs(dt_wl[ok_wl])^2 - 2d*real_part(dt_wl[ok_wl]^2*conj(E[i]))
      Isrc_wl_r[ok_wl]+=((10d)^logS0_r[i]/(2d*!dpi*alpha[i]^2))*$
                        exp(-0.5d*(r2_wl/alpha[i]^2)^(n[i]/2d))
   endif
   if nok_ul gt 0 then begin
      r2_ul=(1d + abs(E[i]^2))*abs(dt_ul[ok_ul])^2 - 2d*real_part(dt_ul[ok_ul]^2*conj(E[i]))
      Isrc_ul_r[ok_ul]+=((10d)^logS0_r[i]/(2d*!dpi*alpha[i]^2))*$
                        exp(-0.5d*(r2_ul/alpha[i]^2)^(n[i]/2d))
   endif
endfor
forloop_status,0,0,wn,/delete
Isrc_ul_b=Isrc_ul_r*(10d)^(-0.4d*bkg_color)
Isrc_wl_b=Isrc_wl_r*(10d)^(-0.4d*bkg_color)

bkgcat[*,0]=real_part(theta_gal)
bkgcat[*,1]=imaginary(theta_gal)
bkgcat[*,2]=real_part(beta_gal)
bkgcat[*,3]=imaginary(beta_gal)
bkgcat[*,4]=logS0_b
bkgcat[*,5]=logS0_r
bkgcat[*,6]=alpha
bkgcat[*,7]=real_part(E)
bkgcat[*,8]=imaginary(E)
bkgcat[*,9]=n


print,'Saving images and data'
save_data,memcat,tag+'_clusmem.dat',comment=memcom
save_data,bkgcat,tag+'_bkggals.dat',comment=bkgcom


; Member Galaxies
mwrfits,Imem_r,tag+'_Imem_r.fits'
mwrfits,Imem_b,tag+'_Imem_b.fits'

; Background Galaxies
mwrfits,Isrc_ul_r,tag+'_Isrc_ul_r.fits'
mwrfits,Isrc_ul_b,tag+'_Isrc_ul_b.fits'
mwrfits,Isrc_wl_r,tag+'_Isrc_wl_r.fits'
mwrfits,Isrc_wl_b,tag+'_Isrc_wl_b.fits'

; Deflection
mwrfits,real_part(delta),tag+'_delta1.fits'
mwrfits,imaginary(delta),tag+'_delta2.fits'

; Observed
psf=mk_gaussian_psf(psffwhm)
Iobs_r=convolve(Isrc_wl_r+Imem_r,psf)+bkg
Iobs_b=convolve(Isrc_wl_b+Imem_b,psf)+bkg
mwrfits,Iobs_r,tag+'_Iobs_nn_r.fits'
mwrfits,Iobs_b,tag+'_Iobs_nn_b.fits'

for i=0,n_elements(img_noise)-1 do begin
   mwrfits,Iobs_r+randomn(seed,imsize)*img_noise[i],tag+'_Iobs_wn'+$
           strcompress(string(i+1),/remove_all)+'_r.fits'
   mwrfits,Iobs_b+randomn(seed,imsize)*img_noise[i],tag+'_Iobs_wn'+$
           strcompress(string(i+1),/remove_all)+'_b.fits'
endfor

end

