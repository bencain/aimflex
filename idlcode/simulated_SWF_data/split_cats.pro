pro split_cats, name, seed=seed, show=show, bcg_ra=bcg_ra, bcg_dec=bcg_dec, $
                wlnoise=wlnoise, flnoise=flnoise, slnoise=slnoise, noise=noise

  if keyword_set(noise) then begin
     wlnoise=1
     slnoise=1
     flnoise=1
  endif

; This procedure splits a simulated image catalog between strong,
; weak, and flexion lensing.  

  sl_cname=name+'_sl.cat'
  wl_cname=name+'_wl.cat'
  fl_cname=name+'_fl.cat'

  sl_dname=name+'_sl.dat'
  wl_dname=name+'_wl.dat'
  fl_dname=name+'_fl.dat'

  cat = read_data(name+'.dat')
  nobj=n_elements(cat[*,0])

  realformat='(E15.8)'          ; need to have it as Ex.y and x >= y+7
  int_format='(I10)'            ; Lots of room

; CAT has these data in the second dimension:
; 0-1   RA,Dec
; 2-3   theta
; 4     m_obs
; 5     r_obs
; 6     z
; 7-8   beta
; 9-10  unlensed RA,Dec
; 11    psi
; 12-13 alpha
; 14-16 kappa, gamma
; 17-20 F, G
; 21    N_family
; 22    N_images
; 23    i_image

  if not keyword_set(bcg_ra)  then bcg_ra= median(cat[*,0])
  if not keyword_set(bcg_dec) then bcg_dec=median(cat[*,1])


; Calculate the reduced shear and flexion
  g=dcomplex(cat[*,15],cat[*,16])/(1d - cat[*,14])
  flex1=dcomplex(cat[*,17],cat[*,18])/(1d - cat[*,14])
  flex3=dcomplex(cat[*,19],cat[*,20])/(1d - cat[*,14])


;;;;;;;;;;
; Deal with noise

; For SL
  sig_theta=0.03d/60d           ; In arcmin

; For WL
  sig_e=0.1d
  
; For FL
  sig_flex1=0.01d               ; In arcsec^-1
  sig_flex3=0.01d
  

  if not keyword_set(slnoise) then dtheta=dblarr(nobj,2) else dtheta=sig_theta*randomn(seed,nobj,2,/double)

  if not keyword_set(wlnoise) then e_src=g*0d else begin
     e_src=sig_e*dcomplex(randomn(seed,nobj,/double),randomn(seed,nobj,/double))
     gt1=where(abs(e_src) gt 1d,ngt1)
     if ngt1 gt 0 then e_src[gt1]=1d/conj(e_src[gt1])
  endelse


;     emag = (-sig_e^2)*alog(!dpi*sig_e^2*(1d - exp(-1d/sig_e^2)) * randomn(seed,nobj,/double)) ; From Schneider 1996
;     phase=2d*!dpi*randomu(seed,nobj,/double)
;     big=where(emag gt 0.9d,nbig)
;     if nbig gt 0 then emag[big]=0.9d
;     e_src=dcomplex(emag*cos(phase),emag*sin(phase))


  if not keyword_set(flnoise) then begin
     intflex1=flex1*0d
     intflex3=flex3*0d
  endif else begin
     intflex1 = sig_flex1*dcomplex(randomn(seed,nobj,/double),randomn(seed,nobj,/double))
     intflex3 = sig_flex3*dcomplex(randomn(seed,nobj,/double),randomn(seed,nobj,/double))
  endelse
  
  devsq_sl=total(dtheta^2)
  devsq_wl=total(abs(e_src)^2)
  devsq_f1=total(abs(intflex1)^2)
  devsq_f3=total(abs(intflex3)^2)
  
  openw,notesunit,name+'_notes.txt',/get_lun,width=18*10,/append
  printf,notesunit,'dev_sl^2 =',devsq_sl
  printf,notesunit,'dev_wl^2 =',devsq_wl
  printf,notesunit,'dev_f1^2 =',devsq_f1
  printf,notesunit,'dev_f3^2 =',devsq_f3
  close,notesunit
  free_lun,notesunit

;;;;;;;;;;

; Create noisy positions - right now this only goes into SL through RADec_obs
  theta_obs=cat[*,2:3] + dtheta ; in arcmin
  RADec_obs=cat[*,0:1]*0d       ; in degrees
  RADec_obs[*,0]=bcg_ra  - (cat[*,2]+dtheta[*,0])/(60d*cos((!dpi/180d)*bcg_dec)) 
  RADec_obs[*,1]=bcg_dec + (cat[*,3]+dtheta[*,1])/60d

; Convert source ellipiticity and shear into observed ellipticity
  gt1=where(abs(g) gt 1d,ngt1,complement=le1,ncomplement=nle1)
  e_img=e_src*0d
  if nle1 gt 0 then e_img[le1]=(e_src[le1] + g[le1])/(dcomplex(1d,0d) + conj(g[le1])*e_src[le1])
  if ngt1 gt 0 then e_img[gt1]=(dcomplex(1d,0d) + conj(e_src[gt1])*g[gt1])/conj(e_src[gt1] + g[gt1])
  
; Convert flexion and flexion noise into observed flexion
  psi1=flex1 + intflex1
  psi3=flex3 + intflex3

; Split up the samples
  flexmax1=0.5d                 ; When is 1-flexion too big to be trusted?
  flexmin1=0d                   ; When is 1-flexion too small to be trusted?

  flexmax3=0.5d                 ; When is 3-flexion too big to be trusted?
  flexmin3=0d                   ; When is 3-flexion too small to be trusted?

  shearlim=0.9d                 ; When is shear too big to be trusted?
  mulim=0.1d                    ; When is a SL source too demagnified?

; Multiple vs single images
  sing_img=where(cat[*,22] eq 1d,n_sing,complement=mult_img,ncomplement=n_mult)

; Now we need to get the SL objects.  "Quad" images have 4 or more
; locations...we'll take the objects with the maximum number of images.
;  slobj=where((cat[*,22] eq max(cat[*,22])) and (cat[*,22] gt 1),n_sl,complement=wlobj,ncomplement=n_wl)
  slobj=where(cat[*,22] gt 1,n_sl,complement=wlobj,ncomplement=n_wl)

  if n_wl gt 0 then begin
     wl_mu=1d/((1d - cat[wlobj,14])^2 - total(cat[wlobj,15:16]^2,2))
     ok=where(abs(wl_mu) gt mulim,nok)
     if nok gt 0 then begin
        wlobj=wlobj[ok]
        n_wl=nok
     endif else begin
        wlobj=-1
        n_wl=0
     endelse
  endif
     


; Get flexion galaxies
;  flobj=where((abs(psi1[wlobj]) lt flexmax) and (abs(psi3[wlobj]) lt flexmax) and $
;              (abs(psi1[wlobj]) gt flexmin) and (abs(psi3[wlobj]) gt flexmin),n_fl)
;  if n_fl gt 0 then flobj=wlobj[flobj]
  flobj=where((abs(psi1) lt flexmax1) and (abs(psi3) lt flexmax3) and $
              (abs(psi1) gt flexmin1) and (abs(psi3) gt flexmin3),n_fl)

  if n_fl gt 0 then begin
     fl_mu=1d/((1d - cat[flobj,14])^2 - total(cat[flobj,15:16]^2,2))
     ok=where(abs(fl_mu) gt mulim,nok)
     if nok gt 0 then begin
        flobj=flobj[ok]
        n_fl=nok
     endif else begin
        flobj=-1
        n_fl=0
     endelse
  endif



; WL and flexion objects below the limits
;  if n_sing gt 0 then begin
;     wlobj=where(abs(e_img[sing_img]) lt shearlim,n_wl)
;     if n_wl gt 0 then wlobj=sing_img[wlobj]
;
;     flobj=where((abs(psi1[sing_img]) lt flexmax) and (abs(psi3[sing_img]) lt flexmax) and $
;                 (abs(psi1[sing_img]) gt flexmin) and (abs(psi3[sing_img]) gt flexmin),n_fl)
;     if n_fl gt 0 then flobj=sing_img[flobj]
;  endif else begin
;     wlobj=-1
;     flobj=-1
;     n_wl=0
;     n_fl=0
;  endelse

  print,'Weak   ',n_wl
  print,'Strong ',n_sl
  print,'Flexion',n_fl

  if (n_wl gt 0) then begin
; Save the WL data in the SWunited format
     wlout=dblarr(n_wl,6)
     wlout[*,0]=cat[wlobj,0]
     wlout[*,1]=cat[wlobj,1]
     wlout[*,2]=real_part(e_img[wlobj])
     wlout[*,3]=imaginary(e_img[wlobj])
     wlout[*,4]=1d
     wlout[*,5]=cat[wlobj,6]

     if keyword_set(show) then begin
        showcat=wlout

        med_ra=median(cat[*,0])
        med_dec=median(cat[*,1])
; Convert to relative arcsec coords
        showcat[*,0]=-cos(!dpi*med_dec/180d)*3600d*(showcat[*,0]-med_ra)
        showcat[*,1]=3600d*(showcat[*,1]-med_dec)

        plot_field,showcat,wnum2,spin=2
     endif

   
     openw,wunit,wl_cname,/get_lun,width=18*6

     printf,wunit,'# WL Data Catalog'
     printf,wunit,'# Coord type, then data'
     printf,wunit,'# RA'
     printf,wunit,'# Dec'
     printf,wunit,'# e1'
     printf,wunit,'# e2'
     printf,wunit,'# weight'
     printf,wunit,'# redshift'
     printf,wunit,'wcs'
     for i=0,n_wl-1 do printf,wunit,$
                              strcompress(string(wlout[i,0],format=realformat),/remove_all)+'  ',$
                              strcompress(string(wlout[i,1],format=realformat),/remove_all)+'  ',$
                              strcompress(string(wlout[i,2],format=realformat),/remove_all)+'  ',$
                              strcompress(string(wlout[i,3],format=realformat),/remove_all)+'  ',$
                              strcompress(string(wlout[i,4],format=realformat),/remove_all)+'  ',$
                              strcompress(string(wlout[i,5],format=realformat),/remove_all)
  
     close,wunit
     free_lun,wunit
  endif
  
  if (n_sl gt 0) then begin
; Save the SL data in the SWunited format
     sl_cat=cat[slobj,*]
     sl_e=e_img[slobj]
     sl_theta=RADec_obs[slobj,*]
     sl_mu=1d/((1d - sl_cat[*,14])^2 - total(sl_cat[*,15:16]^2,2))

; Count the number of systems and the number of images in each system.

     nsys=n_elements(uniq(sl_cat[*,21]))
     nimg=sl_cat[uniq(sl_cat[*,21]),22]
     sys0=total(nimg,/cumulative)-nimg ; The indices in sl_cat of the first image in each system

     order=reverse(sort(sl_cat[sys0,22])) ; Sort them by increasing # of imgs
     sys0=sys0[order]
     nimg=nimg[order]

     openw,wunit,sl_cname,/get_lun,width=18*10

     printf,wunit,'# SL Data Catalog'
     printf,wunit,'# Coord Type'
     printf,wunit,'# Number of SL systems (pimages.nsystems), (pimages.ncur), (pimages.nrec)'
     printf,wunit,'wcs'
     printf,wunit,$
            strcompress(string(nsys,format=int_format),/remove_all)+'  ',$
            strcompress(string(0,format=int_format),/remove_all)+'  ',$
            strcompress(string(0,format=int_format),/remove_all)
     
     printf,wunit,'# Image family, number of images, use for reconstr, source redshift, redshift error'
     printf,wunit,'# Then, for each image in the system: x, y, dx, dy, e1, e2, de1, de2, flux, flux error'
     
     z_sys=sl_cat[sys0,6]

     for i=0,nsys-1 do begin
        ims=sys0[i]+lindgen(nimg[i]) ; Indices to the images in this family
        imuse=where(abs(sl_mu[ims]) gt mulim,nimuse)

        usesys=1

        if nimuse gt 0 then begin
           printf,wunit,'# Image Family '+strcompress(string(sl_cat[sys0[i],21],format=int_format),/remove_all)
           printf,wunit,$
                  strcompress(string(i,format=int_format),/remove_all)+'  ',$        ; System index
                  strcompress(string(nimuse,format=int_format),/remove_all)+'  ',$   ; Number of images
                  strcompress(string(usesys,format=int_format),/remove_all)+'  ',$   ; Whether to use
                  strcompress(string(z_sys[i],format=realformat),/remove_all)+'  ',$ ; System redshift
                  strcompress(string(0.01d,format=realformat),/remove_all)           ; Redshift error

           for j=0,nimuse-1 do begin
              m=ims[imuse[j]]
              printf,wunit,$
                     strcompress(string(sl_theta[m,0],format=realformat),/remove_all)+'  ',$ ; Position
                     strcompress(string(sl_theta[m,1],format=realformat),/remove_all)+'  ',$
                     strcompress(string(sig_theta/60d,format=realformat),/remove_all)+'  ',$      ; Position errors
                     strcompress(string(sig_theta/60d,format=realformat),/remove_all)+'  ',$      ; (in decimal degrees)
                     strcompress(string(real_part(sl_e[m]),format=realformat),/remove_all)+'  ',$ ; Ellipticity
                     strcompress(string(imaginary(sl_e[m]),format=realformat),/remove_all)+'  ',$
                     strcompress(string(0.01d,format=realformat),/remove_all)+'  ',$ ; Ellipticity errors
                     strcompress(string(0.01d,format=realformat),/remove_all)+'  ',$ 
                     strcompress(string(sl_mu[m],format=realformat),/remove_all)+'  ',$ ; Flux - isn't used, so put in mu for ref.
                     strcompress(string(0.1d,format=realformat),/remove_all)            ; Flux error
           endfor
        endif
     endfor
     
     
     
     close,wunit
     free_lun,wunit
  endif


  if (n_fl gt 0) then begin
; Save the flexion data into the SWunited format
     openw,wunit,fl_cname,/get_lun,width=18*10

     flout=dblarr(n_fl,10)
     flout[*,0]=cat[flobj,0]
     flout[*,1]=cat[flobj,1]
     flout[*,2]=real_part(e_img[flobj])
     flout[*,3]=imaginary(e_img[flobj])
     flout[*,4]=real_part(psi1[flobj])
     flout[*,5]=imaginary(psi1[flobj])
     flout[*,6]=real_part(psi3[flobj])
     flout[*,7]=imaginary(psi3[flobj])
     flout[*,8]=1d
     flout[*,9]=cat[flobj,6]
     

     if keyword_set(show) then begin
        showcat=flout

        med_ra=median(cat[*,0])
        med_dec=median(cat[*,1])
; Convert to relative arcsec coords
        showcat[*,0]=-cos(!dpi*med_dec/180d)*3600d*(showcat[*,0]-med_ra)
        showcat[*,1]=3600d*(showcat[*,1]-med_dec)

        plot_field,showcat,wnum1,spin=1,field_cols=[4,5],/log
        plot_field,showcat,wnum3,spin=3,field_cols=[6,7],/log
     endif


     printf,wunit,'# Flexion Data Catalog'
     printf,wunit,'# Coord Type, then data'
     printf,wunit,'# RA'
     printf,wunit,'# Dec'
     printf,wunit,'# e1'
     printf,wunit,'# e2'
     printf,wunit,'# Psi11'
     printf,wunit,'# Psi12'
     printf,wunit,'# Psi31'
     printf,wunit,'# Psi32'
     printf,wunit,'# weight'
     printf,wunit,'# redshift'
     printf,wunit,'wcs'
     for i=0,n_fl-1 do printf,wunit,$
                              strcompress(string(flout[i,0],format=realformat),/remove_all)+'  ',$
                              strcompress(string(flout[i,1],format=realformat),/remove_all)+'  ',$
                              strcompress(string(flout[i,2],format=realformat),/remove_all)+'  ',$
                              strcompress(string(flout[i,3],format=realformat),/remove_all)+'  ',$
                              strcompress(string(flout[i,4],format=realformat),/remove_all)+'  ',$
                              strcompress(string(flout[i,5],format=realformat),/remove_all)+'  ',$
                              strcompress(string(flout[i,6],format=realformat),/remove_all)+'  ',$
                              strcompress(string(flout[i,7],format=realformat),/remove_all)+'  ',$
                              strcompress(string(flout[i,8],format=realformat),/remove_all)+'  ',$
                              strcompress(string(flout[i,9],format=realformat),/remove_all)
     close,wunit
     free_lun,wunit
  endif


; Save the true and observed values of everything
; CAT has these data in the second dimension:
; 0-1   RA,Dec
; 2-3   theta
; 4     m_obs
; 5     r_obs
; 6     z
; 7-8   unlensed RA,Dec
; 9-10  beta
; 11    psi
; 12-13 alpha
; 14-16 kappa, gamma
; 17-20 F, G
; 21    N_family
; 22    N_images
; 23    i_image

; TRUECAT has these data in the second dimension:
; 0-1   RA,Dec
; 2-3   theta
; Galaxy properties
; 4     m_obs 
; 5     r_obs
; 6     z
; 7-8   ellipticity
; Unlensed position
; 9-10  unlensed RA,Dec
; 11-12 beta
; True lensing fields
; 13    psi
; 14-15 alpha
; 16-18 kappa, gamma
; 19-22 F, G
; Lensing Estimators
; 23-24 Observed image position
; 25-26 Observed ellipticity
; 27-28 Observed 1-flexion
; 29-30 Observed 3-flexion
; Image families
; 31    N_family
; 32    N_images
; 33    i_image


  truecat=dblarr(nobj,34)
  truecat[*,0:6]=cat[*,0:6]      ; Positions and intrinsic properties
  truecat[*,7]=real_part(e_src)  ; Unlensed ellipticity
  truecat[*,8]=imaginary(e_src)  ;
  truecat[*,9:12]=cat[*,7:10]    ; Unlensed RA,Dec,beta1,2
  truecat[*,13:22]=cat[*,11:20]  ; True lensing fields
  truecat[*,23]=RADec_obs[*,0]   ; SL position estimator
  truecat[*,24]=RADec_obs[*,1]   ;
  truecat[*,25]=real_part(e_img) ; Shear estimator
  truecat[*,26]=imaginary(e_img) ;
  truecat[*,27]=real_part(psi1)  ; 1-flexion estimator
  truecat[*,28]=imaginary(psi1)  ;
  truecat[*,29]=real_part(psi3)  ; 3-flexion estimator
  truecat[*,30]=imaginary(psi3)  ;
  truecat[*,31]=cat[*,21]        ; N_family
  truecat[*,32]=cat[*,22]        ; N_images
  truecat[*,33]=cat[*,23]        ; i_image

  comm='RA,Dec; theta1,2; M,R,z; e_src1,2; unlensedRA,Dec; beta1,2; psi; alpha1,2; kappa; gamma1,2; F1,2; G1,2; '+$
       'RA_obs,Dec_obs; e_img1,2; psi11,2; psi31,2; N_fam,N_imgs,i_img'

  save_data,truecat,name+'_true.dat',comment=comm

  save_data,truecat[wlobj,*],wl_dname,comment=comm
  save_data,truecat[slobj,*],sl_dname,comment=comm
  save_data,truecat[flobj,*],fl_dname,comment=comm

end

