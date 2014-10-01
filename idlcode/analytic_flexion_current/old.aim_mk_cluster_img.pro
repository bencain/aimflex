FUNCTION AIM_MK_CLUSTER_IMG, IMSIZE, NGAL, BKG, TAG,$
                            IMG_NOISE=IMG_NOISE,$
                            SERSIC=SERSIC, MOFFAT=MOFFAT,$
                            PSEUDOGAUSS=PSEUDOGAUSS, $
                            NFW_SCALE=NFW_SCALE,$
                            PSFFWHM=PSFFWHM

starttime=systime()
if not keyword_set(img_noise) then img_noise=1d-6

; Full image base
image=dblarr(imsize,imsize)+bkg

; Galaxy positions
nside=ceil(sqrt(ngal))
gap=double(imsize)/double(nside+1)
galx=gap*(dindgen(ngal) mod nside)+gap
galy=gap*(double(lindgen(ngal)/long(nside)))+gap
pos=dblarr(ngal,2)
pos[*,0]=galx
pos[*,1]=galy

; Copy the parameter draws from AIM_MCMC_EXPLORE
if not keyword_set(psffwhm) then psffwhm=1.8d
if (keyword_set(sersic) or keyword_set(moffat)) then npar=13 else $
   if keyword_set(pseudogauss) then npar=14 else npar=12

trueps=dblarr(ngal,npar+3) ; An index number, field position, plus true parameters

; Filenames for the output
truefile= tag+'_true.dat'

; Output space
tp=dblarr(npar)



cutrad_scale=25d

; Set up the parameter names
tpnames=strarr(npar)
tpnames[0:5]=['logN0','X_c','Y_c','Alpha','e_+','e_x']
tpnames[npar-6:npar-1]=['g_1','g_2','psi1_1','psi1_2','psi3_1','psi3_2']
if keyword_set(sersic) then tpnames[6]='n'
if keyword_set(moffat) then tpnames[6]='b'
if keyword_set(pseudogauss) then tpnames[6:7]=['k1','k2']

bkg=0.07d                       ; Doing things in flux units

; Get the true parameter ranges.
pmins=dblarr(npar)
pmaxs=dblarr(npar)

pmins[0:5]=[-1d,$                ;logN0
            -5d,-5d,$           ;X_c, Y_c
            3d,$                ;alpha
            -0.8d,-0.8d]        ;Eplus, Ecross

pmins[npar-6:npar-1]=[-0.8d,-0.8d,$           ;g1, g2
                      -2d-2,-2d-2,$           ;psi11,psi12
                      -2d-2,-2d-2]            ;psi31,psi32

pmaxs[0:5]=[1d,$                ;logN0
            5d,5d,$             ;X_c, Y_c
            15d,$               ;alpha
            0.8d,0.8d]          ;Eplus, Ecross

pmaxs[npar-6:npar-1]=[0.8d,0.8d,$           ;g1, g2
                      2d-2,2d-2,$           ;psi11,psi12
                      2d-2,2d-2]            ;psi31,psi32

; So that the peak is at least NSIG sigma above the background level (the
; max bit is in there so that the min isn't TOO small).
nsig=1.5d
if img_noise gt 0d then $
   pmins[0]=max([alog10(2d*!dpi*pmaxs[3]^2*bkg*img_noise*nsig),pmins[0]])


if keyword_set(sersic) then begin
; The Sersic index will run from 0.25 to 4
   pmins[6]=0.25d
   pmaxs[6]=4d
endif else if keyword_set(pseudogauss) then begin
; The ki parameters will run from -5 to 5
   pmins[6:7]=[0d,0d]
   pmaxs[6:7]=[5d,5d]
endif else if keyword_set(moffat) then begin
; The Moffat slope will run from close to 1 to 5 or so
   pmins[6]=1.2d
   pmaxs[6]=7d
endif

if keyword_set(fixcenter) then begin
   pmins[1:2]=0d
   pmaxs[1:2]=0d
endif


; Make some images and put them into the big image
for i=0,ngal-1 do begin

   forloop_status,i,ngal,flnum,label='Making galaxy images...'
 
   tp=aim_parameter_draw(pmins,pmaxs,seed,ctr_lim=pmaxs[1],e_lim=pmaxs[4],$
                        sersic=keyword_set(sersic),$
                        pseudogauss=keyword_set(pseudogauss),$
                        moffat=keyword_set(moffat))

   if not keyword_set(nfw_scale) then nfw_scale=1d
   z_lens=0.183d
   r200=2.16d                   ; in h_70^-1 Mpc (Limousin et al. 2007)
   c200=7.6d                    ; also Limousin et al. 2007
;      r200=2.4d                                 ; in h_70^-1 Mpc (Coe et al. 2010)
;      c200=9.2d                                 ; also Coe et al. 2010
   d_ang=cosmo_dist(0.183)                          ; Distance to A1689 in h_70^-1 Mpc
   theta200=(r200/d_ang)*(180d/!dpi)*3600d/0.05d    ; in HST pixels

   theta_s=theta200/c200
   lpars=nfw_lensmodel(pos[i,*],ctr=[imsize/2d,imsize/2d],$
                       theta_s=theta_s,z=z_lens,c200=c200,force_scale=nfw_scale)
   tp[npar-6:npar-1]=lpars

   trueps[i,*] =[i,pos[i,0],pos[i,1],tp]

   cutrad=double(ceil(cutrad_scale*tp[3]))
   win=mk_window(ceil((2d)*cutrad))
   psf=mk_gaussian_psf(psffwhm)

   dataimg=aim_mk_image(tp,0d,win,psf,$
                           sersic=keyword_set(sersic),$
                           moffat=keyword_set(moffat),$
                           pseudogauss=keyword_set(pseudogauss))
   

   image=insert_stamp(image,dataimg,pos[i,*],/sum)
;   mwrfits,dataimg,tag+strcompress(string(i),/remove_all)+'.fits'

endfor

image*=(1d + img_noise*randomn(seed,imsize,imsize))

trueparlist='logN0, Xc, Yc, alpha, Eplus, Ecross, '
if keyword_set(sersic_true) then trueparlist+='n, '
if keyword_set(pseudogauss_true) then trueparlist+='k1, k2, '
if keyword_set(moffat_true) then trueparlist+='b, '
trueparlist+='g_1, g_2, psi1_1, psi1_2, psi3_1, psi3_2'

tpcom='True parameter values: N, X_field, Y_field, '+trueparlist
save_data,trueps,truefile,comment=tpcom
forloop_status,0,0,flnum,/delete

print,'Started:  '+starttime
print,'Finished: '+systime()

return,image

end
