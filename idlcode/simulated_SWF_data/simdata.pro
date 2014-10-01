pro simdata, name, density, fieldsize, nrealization, $
             z_lens=z_lens, seed=seed, $
             n_slsys=n_slsys, sl_fieldsize=sl_fieldsize, $
             header=header, bcg_ra=bcg_ra,bcg_dec=bcg_dec,$
             lens0=lens0, lens1=lens1, lens2=lens2, lens3=lens3, lens4=lens4, lens5=lens5,$
             noiseoff=noiseoff

  if not keyword_set(n_slsys) then n_slsys=3L
  if n_slsys lt 3 then n_slsys=3L
  if keyword_set(noiseoff) then nflag=0d else nflag=1d

  if not keyword_set(sl_fieldsize) then sl_fieldsize=0.2d
  sl_fieldsize=fieldsize/4d


; Allow for a fits file header to be copied
  if not keyword_set(header) then header='/Users/bcain/idl/pro/simulated_SWF_data/1e0657R_avg_wcs_hdr.fits'
  hdr=headfits(header)

  naxis1=sxpar(hdr,'NAXIS1')
  naxis2=sxpar(hdr,'NAXIS2')
  
  c1=naxis1/2d
  c2=naxis2/2d

  realformat='(E15.8)'          ; need to have it as Ex.y and x >= y+7
  int_format='(I10)'            ; Lots of room

  ; Set the RA/Dec of the BCG to the field center if not otherwise specified
  if not (keyword_set(bcg_ra) and keyword_set(bcg_dec)) then $
     xyad,hdr,c1,c2,bcg_ra,bcg_dec,/celestial,/print

; How many objects do we want in the final catalog
  nsrc=ceil(density*fieldsize^2,/l64)

; Lens model details
  if not keyword_set(z_lens) then z_lens=0.2d
  default=[0d,0d,$              ; Center
           0d,0d,1d,$           ; sigma, k0, theta_c
           0d,0d,1d]            ; eps, phi, DLInf/DInf

  if n_elements(lens0) lt 8 then begin
     lens0=default
     lens0[2]=1d
     lens0[4]=4d-2
     lens0[5]=0.0d
     lens0[7]=redshift_to_weight(10000d,z_lens)
  endif

  haloes=lindgen(6)
  lens=lens0

  if not keyword_set(lens1) then begin
     lens1=default
     haloes=set_difference(haloes,1L)
  endif else lens=[lens,lens1]
  if not keyword_set(lens2) then begin
     lens2=default
     haloes=set_difference(haloes,2L)
  endif else lens=[lens,lens2]
  if not keyword_set(lens3) then begin
     lens3=default
     haloes=set_difference(haloes,3L)
  endif else lens=[lens,lens3]
  if not keyword_set(lens4) then begin
     lens4=default
     haloes=set_difference(haloes,4L)
  endif else lens=[lens,lens4]
  if not keyword_set(lens5) then begin
     lens5=default
     haloes=set_difference(haloes,5L)
  endif else lens=[lens,lens5]
  
  nhalo=n_elements(haloes)
  for i=0,nhalo-1 do print,"Halo at ",lens[8*i],lens[8*i+1]," with sigma = ",lens[8*i+2]

; Fitting tolerance in the source plane
  errs=1d-12                    ; in arcmin
  erri=1d/60d                   ; in arcmin

; For splitting up the samples
  flexmax1=1d                   ; When is 1-flexion too big to be trusted?
  flexmin1=0.00d               ; When is 1-flexion too small to be trusted?
  
  flexmax3=1d                   ; When is 3-flexion too big to be trusted?
  flexmin3=0.00d               ; When is 3-flexion too small to be trusted?
  
  shearlim=0.95d                ; When is shear too big to be trusted?

  mulim=0.9d                   ; When is a SL source too demagnified?
  if nhalo gt 1 then nimmin=4 else nimmin=2



; The catalog output will have the following entries
; 0-1   RA, Dec
; 2-3   theta1, theta2
; 4-5   z, Z_wt
; 6-7   RA_src, Dec_src
; 8-9   beta1, beta2
; 10    psi
; 11-12 alpha1, alpha2
; 13    kappa
; 14-15 gamma1, gamma2
; 16-17 fflex1, fflex2
; 18-19 gflex1, gflex2
; 20-22 N_family, N_images, i_image

  cat=dblarr(nsrc,23)

; Now get the SL systems.  We'll draw source positions until
; the number of multiple image systems is obtained.
  
  nok=0
  count=0
  family=0

; Set up an output plot for the SL images

  print, 'Finding SL images'
  while nok lt n_slsys do begin
     
     bok=0
     while not bok do begin
        beta = sl_fieldsize*(randomu(seed,2,/double) - 0.5d) ; Choose a source position
        if sqrt(total(beta^2)) lt sl_fieldsize then bok=1
     endwhile
     
; Pick which halo to center around
;     which=haloes[floor(randomu(seed)*nhalo)]
;     print,which
;     if which eq 1 then beta+=lens1[0:1]
;     if which eq 2 then beta+=lens2[0:1]
;     if which eq 3 then beta+=lens3[0:1]
;     if which eq 4 then beta+=lens4[0:1]
;     if which eq 5 then beta+=lens5[0:1]

     zedok=0
     z0=2d/3d
     k=3L
     while not zedok do begin
        zed = z0*randomu(seed,/double,gamma=k) ; Get the tentative redshift
        if zed gt z_lens then zedok=1          ; Only accept redshifts beyond the 
     endwhile
     
     Zwt = redshift_to_weight(zed, z_lens)                ; Get the weight
     thetas = find_mult_images_triangles(beta, errs, zwt, lens, $
                                         span=2d*fieldsize, ngrid=20L, $
                                         nimages=nimages, resolution=erri) ; Find the images

     lens_fields=my_lens(thetas,lens,z_weight=Zwt)
     mu=1d/((1d - lens_fields[*,3])^2-total(lens_fields[*,4:5]^2,2))

     obs=where(abs(mu) gt mulim,nobs)


; If it's a multiple image system, then we'll record it
     if nobs ge nimmin then begin

        plot,[-0.5,0.5,0.5,-0.5,-0.5]*fieldsize,[-0.5,-0.5,0.5,0.5,-0.5]*fieldsize,/isotropic
        oplot,[lens[0]],[lens[1]],psym=4,symsize=3,color=fsc_color('red'),thick=3
        for i=1,nhalo-1 do oplot,[lens[8*i]],[lens[8*i+1]],psym=4,symsize=3,color=fsc_color('blue'),thick=2
        oplot,thetas[obs,0],thetas[obs,1],/psym
        for i=0,nobs-1 do oplot,[thetas[obs[i],0],beta[0]],[thetas[obs[i],1],beta[1]],color=fsc_color('green')
        
        print,'beta = ',beta[0],beta[1]
        read,'Is this system ok? (1=y, 0=no)',sysok

        if sysok eq 1 then begin
           cat[count:count+nobs-1, 0]=bcg_ra - thetas[obs,0]/(60d*cos((!dpi/180d)*bcg_dec)) ; RA - observed
           cat[count:count+nobs-1, 1]=bcg_dec + thetas[obs,1]/60d                           ; Dec -observed
           cat[count:count+nobs-1, 2]=thetas[obs,0]                                         ; Theta1
           cat[count:count+nobs-1, 3]=thetas[obs,1]                                         ; Theta2
           cat[count:count+nobs-1, 4]=zed                                                   ; redshift
           cat[count:count+nobs-1, 5]=zwt                                                   ; redshift weight
           cat[count:count+nobs-1, 6]=bcg_ra - beta[0]/(60d*cos((!dpi/180d)*bcg_dec))       ; RA - unlensed
           cat[count:count+nobs-1, 7]=bcg_dec+ beta[1]/60d                                  ; Dec -unlensed
           cat[count:count+nobs-1, 8]=beta[0]                                               ; beta1
           cat[count:count+nobs-1, 9]=beta[1]                                               ; beta2
           cat[count:count+nobs-1,10]=lens_fields[obs,0]                                    ; psi
           cat[count:count+nobs-1,11]=lens_fields[obs,1]                                    ; alpha1
           cat[count:count+nobs-1,12]=lens_fields[obs,2]                                    ; alpha2
           cat[count:count+nobs-1,13]=lens_fields[obs,3]                                    ; kappa
           cat[count:count+nobs-1,14]=lens_fields[obs,4]                                    ; gamma1
           cat[count:count+nobs-1,15]=lens_fields[obs,5]                                    ; gamma2
           cat[count:count+nobs-1,16]=lens_fields[obs,6]/60d                                ; fflex1   Note that the flexion is calculated
           cat[count:count+nobs-1,17]=lens_fields[obs,7]/60d                                ; fflex2   in units of 1/arcmin, but we change
           cat[count:count+nobs-1,18]=lens_fields[obs,8]/60d                                ; gflex1   that to 1/arcsec here.
           cat[count:count+nobs-1,19]=lens_fields[obs,9]/60d                                ; gflex2
           cat[count:count+nobs-1,20]=family                                                ; Image family
           cat[count:count+nobs-1,21]=nobs                                                  ; Number of images
           cat[count:count+nobs-1,22]=dindgen(nobs)                                         ; Image index
           
           family++
           count+=nobs
           nok++
           print,nok
        endif
     endif
  endwhile

; Now for the singly lensed systems
  print, 'Finding WL/FL images'
  nsingle=nsrc - count

  thetas=fieldsize*(randomu(seed,nsingle,2,/double) - 0.5d)

  zed = dblarr(nsingle)
  i=0
  z0=2d/3d
  k=3L
  while i lt nsingle do begin
     zee = z0*randomu(seed,/double,gamma=k) ; Get the tentative redshift
     if zee gt z_lens then begin
        zed[i]=zee              ; Only accept redshifts beyond the lens
        i++
     endif
  endwhile

  Zwt = redshift_to_weight(zed, z_lens)       ; Get the weight
  lens_fields=my_lens(thetas,lens,z_weight=Zwt)

  beta=thetas - lens_fields[*,1:2]

; starting with here...
  cat[count:nsrc-1, 0]=bcg_ra - thetas[*,0]/(60d*cos((!dpi/180d)*bcg_dec)) ; RA - observed
  cat[count:nsrc-1, 1]=bcg_dec + thetas[*,1]/60d                           ; Dec -observed
  cat[count:nsrc-1, 2]=thetas[*,0]                                         ; Theta1
  cat[count:nsrc-1, 3]=thetas[*,1]                                         ; Theta2
  cat[count:nsrc-1, 4]=zed                                                 ; redshift
  cat[count:nsrc-1, 5]=zwt                                                 ; redshift weight
  cat[count:nsrc-1, 6]=bcg_ra - beta[*,0]/(60d*cos((!dpi/180d)*bcg_dec))   ; RA - unlensed
  cat[count:nsrc-1, 7]=bcg_dec+ beta[*,1]/60d                              ; Dec -unlensed
  cat[count:nsrc-1, 8]=beta[*,0]                                           ; beta1
  cat[count:nsrc-1, 9]=beta[*,1]                                           ; beta2
  cat[count:nsrc-1,10]=lens_fields[*,0]                                    ; psi
  cat[count:nsrc-1,11]=lens_fields[*,1]                                    ; alpha1
  cat[count:nsrc-1,12]=lens_fields[*,2]                                    ; alpha2
  cat[count:nsrc-1,13]=lens_fields[*,3]                                    ; kappa
  cat[count:nsrc-1,14]=lens_fields[*,4]                                    ; gamma1
  cat[count:nsrc-1,15]=lens_fields[*,5]                                    ; gamma2
  cat[count:nsrc-1,16]=lens_fields[*,6]/60d                                ; fflex1   Note that the flexion is calculated
  cat[count:nsrc-1,17]=lens_fields[*,7]/60d                                ; fflex2   in units of 1/arcmin, but we change
  cat[count:nsrc-1,18]=lens_fields[*,8]/60d                                ; gflex1   that to 1/arcsec here.
  cat[count:nsrc-1,19]=lens_fields[*,9]/60d                                ; gflex2
  cat[count:nsrc-1,20]=dindgen(nsingle)+family                             ; Image family
  cat[count:nsrc-1,21]=1d                                                  ; Number of images
  cat[count:nsrc-1,22]=0d                                                  ; Image index
  
        
; Center the field on the BCG (WCS coords in decimal degrees)
; Keep the left/right orientation correct w.r.t. RA, and account for
; the declination.

  save_data,cat,name+'.dat',comment=' RA, Dec, theta1, theta2, redshift, Z(z), '+$
            'beta_RA, beta_Dec, beta1, beta2, psi, alpha1, alpha2, kappa, gamma1, gamma2, F1, F2, G1, G2, '+$
            'Image_family, N_images, i_image'

  field=[-0.5d,0.5d]*fieldsize

; Save some images of the lensing fields for a source at infinity.
  if keyword_set(header) then begin
                                
; Find the edges of the field in RA/Dec space 
     dec_1=bcg_dec+fieldsize/(2d*60d)
     dec_0=bcg_dec-fieldsize/(2d*60d)

     print,'Find field edges'
     adxy,hdr,bcg_ra,dec_1,x1,y1,/PRINT
     adxy,hdr,bcg_ra,dec_0,x0,y0,/PRINT

     ny=floor(y1-y0,/l64)
     nx=ny

     print,'FITS size'
     print,nx,ny

     sxaddpar,hdr,'NAXIS1',nx
     sxaddpar,hdr,'NAXIS2',ny
     sxaddpar,hdr,'CRPIX1',(nx-1)/2
     sxaddpar,hdr,'CRPIX2',(ny-1)/2
     sxaddpar,hdr,'CRVAL1',bcg_ra
     sxaddpar,hdr,'CRVAL2',bcg_dec

  endif else begin
     nx=floor(fieldsize/(0.1d/60d),/l64)
     ny=nx
  endelse

  thetapix=dblarr(nx*ny,2)
  thetapix[*,0]=fieldsize*((dindgen(nx*ny) mod nx) - (nx-1)/2d)/(nx-1)
  thetapix[*,1]=fieldsize*((lindgen(nx*ny)/long(nx)) - (ny-1)/2d)/(ny-1)
  
  lenspix=my_lens(thetapix,lens)

  image=dblarr(nx,ny)
  
  image[*]=lenspix[*,0]
  mwrfits,image,name+'_psi.fits',hdr

  image[*]=lenspix[*,1]
  mwrfits,image,name+'_alpha1.fits',hdr
  image[*]=lenspix[*,2]
  mwrfits,image,name+'_alpha2.fits',hdr

  image[*]=lenspix[*,3]
  mwrfits,image,name+'_kappa.fits',hdr
  image[*]=lenspix[*,4]
  mwrfits,image,name+'_gamma1.fits',hdr
  image[*]=lenspix[*,5]
  mwrfits,image,name+'_gamma2.fits',hdr

  image[*]=1d/((1d - lenspix[*,3])^2 - lenspix[*,4]^2 - lenspix[*,5]^2)
  hi=where(image[*] gt 100d,nhi)
  lo=where(image[*] lt -100d,nlo)
  if nhi gt 0 then image[hi]=100d
  if nlo gt 0 then image[lo]=-100d
  mwrfits,image,name+'_mu.fits',hdr

  image[*]=lenspix[*,6]/60d
  mwrfits,image,name+'_fflex1.fits',hdr
  image[*]=lenspix[*,7]/60d
  mwrfits,image,name+'_fflex2.fits',hdr
  image[*]=lenspix[*,8]/60d
  mwrfits,image,name+'_gflex1.fits',hdr
  image[*]=lenspix[*,9]/60d
  mwrfits,image,name+'_gflex2.fits',hdr


; Start a notes file
  openw,notesunit,name+'_notes.txt',/get_lun,width=18*10
  printf,notesunit,'header = ',header
  printf,notesunit,'BCG_RA  =',bcg_ra
  printf,notesunit,'BCG_DEC =',bcg_dec
  printf,notesunit,'lens0 = ',lens0
  printf,notesunit,'lens1 = ',lens1
  printf,notesunit,'lens2 = ',lens2
  printf,notesunit,'lens3 = ',lens3
  printf,notesunit,'lens4 = ',lens4
  printf,notesunit,'lens5 = ',lens5
  printf,notesunit,'z_lens = ',z_lens
  printf,notesunit,'nsrc = ',nsrc
  printf,notesunit,'nrealization = ',nrealization


; Define some properties for all the output data catalogs
  g=dcomplex(cat[*,14],cat[*,15])/(1d - cat[*,13])
  flex1=dcomplex(cat[*,16],cat[*,17])/(1d - cat[*,13])
  flex3=dcomplex(cat[*,18],cat[*,19])/(1d - cat[*,13])


; For SL
  sig_theta=3d-1/60d            ; In arcmin, the true errors
  sl_err=0.5d/3600              ; What is the sl error for the reconstruction (in decimal degrees)

; For WL
  sig_e=0.2d
  
; For FL
  sig_flex1=0.01d               ; In arcsec^-1
  sig_flex3=0.01d



; Now we need data catalogs for each of the realizations specified.

  for realization=0,nrealization-1 do begin
     tag=strcompress('_'+string(realization,format='(I03)'),/remove_all)

     sl_cname=name+tag+'_sl.cat'
     wl_cname=name+tag+'_wl.cat'
     fl_cname=name+tag+'_fl.cat'
     
     sl_dname=name+tag+'_sl.dat'
     wl_dname=name+tag+'_wl.dat'
     fl_dname=name+tag+'_fl.dat'
  

; Lens model file
     openw,wunit,name+tag+'_lensmodel',/get_lun,width=18*10
  
     printf,wunit,'# Lens model parameter file, with BCG and other stuff also'
     printf,wunit,'wcs'
  
     printf,wunit,strcompress(string(1,format=int_format),/remove_all) ; 1 lens

     printf,wunit,$
            strcompress(string(2,format=int_format),/remove_all)+'  ',$        ; NIE lens
            strcompress(string(z_lens,format=realformat),/remove_all)+'  ',$   ; Lens redshift
            strcompress(string(bcg_ra,format=realformat),/remove_all)+'  ',$   ; x_ctr position
            strcompress(string(bcg_dec,format=realformat),/remove_all)+'  ',$  ; y_ctr position
            strcompress(string(lens0[2],format=realformat),/remove_all)+'  ',$ ; Velocity dispersion of the NIE in 10^3 km/s
            strcompress(string(lens0[3],format=realformat),/remove_all)+'  ',$ ; constant kappa offset
; We'll pretend that we're ignorant to any ellipticity or substructure.
            strcompress(string(lens0[4],format=realformat),/remove_all)+'  ',$ ; core radius (in arcmin)
            strcompress(string(lens0[5]*0d,format=realformat),/remove_all)+'  ',$ ; ellipticity (a-b)/(a+b)
            strcompress(string(lens0[6]*0d,format=realformat),/remove_all)        ; ellipse position angle

     printf,wunit,'# True model parameters'
     printf,wunit,'# Halo 0: '+$
            strcompress(string(2,format=int_format),/remove_all)+'  ',$     ; NIE lens
            strcompress(string(z_lens,format=realformat),/remove_all)+'  ',$ ; Lens redshift
            strcompress(string(bcg_ra,format=realformat),/remove_all)+'  ',$   ; x_ctr position
            strcompress(string(bcg_dec,format=realformat),/remove_all)+'  ',$  ; y_ctr position
            strcompress(string(lens0[2],format=realformat),/remove_all)+'  ',$ ; Velocity dispersion of the NIE in 10^3 km/s
            strcompress(string(lens0[3],format=realformat),/remove_all)+'  ',$ ; constant kappa offset
            strcompress(string(lens0[4],format=realformat),/remove_all)+'  ',$ ; core radius (in arcmin)
            strcompress(string(lens0[5],format=realformat),/remove_all)+'  ',$ ; ellipticity (a-b)/(a+b)
            strcompress(string(lens0[6],format=realformat),/remove_all)        ; ellipse position angle
     printf,wunit,'# Halo 1: '+$
            strcompress(string(2,format=int_format),/remove_all)+'  ',$     ; NIE lens
            strcompress(string(z_lens,format=realformat),/remove_all)+'  ',$ ; Lens redshift
            strcompress(string(bcg_ra - lens1[0]/(60d*cos((!dpi/180d)*bcg_dec)),format=realformat),/remove_all)+'  ',$ ; x_ctr position
            strcompress(string(bcg_dec+ lens1[1]/60d,format=realformat),/remove_all)+'  ',$ ; y_ctr position
            strcompress(string(lens1[2],format=realformat),/remove_all)+'  ',$ ; Velocity dispersion of the NIE in 10^3 km/s
            strcompress(string(lens1[3],format=realformat),/remove_all)+'  ',$ ; constant kappa offset
            strcompress(string(lens1[4],format=realformat),/remove_all)+'  ',$ ; core radius (in arcmin)
            strcompress(string(lens1[5],format=realformat),/remove_all)+'  ',$ ; ellipticity (a-b)/(a+b)
            strcompress(string(lens1[6],format=realformat),/remove_all)        ; ellipse position angle
     printf,wunit,'# Halo 2: '+$
            strcompress(string(2,format=int_format),/remove_all)+'  ',$     ; NIE lens
            strcompress(string(z_lens,format=realformat),/remove_all)+'  ',$ ; Lens redshift
            strcompress(string(bcg_ra - lens2[0]/(60d*cos((!dpi/180d)*bcg_dec)),format=realformat),/remove_all)+'  ',$ ; x_ctr position
            strcompress(string(bcg_dec+ lens2[1]/60d,format=realformat),/remove_all)+'  ',$ ; y_ctr position
            strcompress(string(lens2[2],format=realformat),/remove_all)+'  ',$ ; Velocity dispersion of the NIE in 10^3 km/s
            strcompress(string(lens2[3],format=realformat),/remove_all)+'  ',$ ; constant kappa offset
            strcompress(string(lens2[4],format=realformat),/remove_all)+'  ',$ ; core radius (in arcmin)
            strcompress(string(lens2[5],format=realformat),/remove_all)+'  ',$ ; ellipticity (a-b)/(a+b)
            strcompress(string(lens2[6],format=realformat),/remove_all)        ; ellipse position angle
     printf,wunit,'# Halo 3: '+$
            strcompress(string(2,format=int_format),/remove_all)+'  ',$     ; NIE lens
            strcompress(string(z_lens,format=realformat),/remove_all)+'  ',$ ; Lens redshift
            strcompress(string(bcg_ra - lens3[0]/(60d*cos((!dpi/180d)*bcg_dec)),format=realformat),/remove_all)+'  ',$ ; x_ctr position
            strcompress(string(bcg_dec+ lens3[1]/60d,format=realformat),/remove_all)+'  ',$ ; y_ctr position
            strcompress(string(lens3[2],format=realformat),/remove_all)+'  ',$ ; Velocity dispersion of the NIE in 10^3 km/s
            strcompress(string(lens3[3],format=realformat),/remove_all)+'  ',$ ; constant kappa offset
            strcompress(string(lens3[4],format=realformat),/remove_all)+'  ',$ ; core radius (in arcmin)
            strcompress(string(lens3[5],format=realformat),/remove_all)+'  ',$ ; ellipticity (a-b)/(a+b)
            strcompress(string(lens3[6],format=realformat),/remove_all)        ; ellipse position angle
     printf,wunit,'# Halo 4: '+$
            strcompress(string(2,format=int_format),/remove_all)+'  ',$     ; NIE lens
            strcompress(string(z_lens,format=realformat),/remove_all)+'  ',$ ; Lens redshift
            strcompress(string(bcg_ra - lens4[0]/(60d*cos((!dpi/180d)*bcg_dec)),format=realformat),/remove_all)+'  ',$ ; x_ctr position
            strcompress(string(bcg_dec+ lens4[1]/60d,format=realformat),/remove_all)+'  ',$ ; y_ctr position
            strcompress(string(lens4[2],format=realformat),/remove_all)+'  ',$ ; Velocity dispersion of the NIE in 10^3 km/s
            strcompress(string(lens4[3],format=realformat),/remove_all)+'  ',$ ; constant kappa offset
            strcompress(string(lens4[4],format=realformat),/remove_all)+'  ',$ ; core radius (in arcmin)
            strcompress(string(lens4[5],format=realformat),/remove_all)+'  ',$ ; ellipticity (a-b)/(a+b)
            strcompress(string(lens4[6],format=realformat),/remove_all)        ; ellipse position angle
     printf,wunit,'# Halo 5: '+$
            strcompress(string(2,format=int_format),/remove_all)+'  ',$     ; NIE lens
            strcompress(string(z_lens,format=realformat),/remove_all)+'  ',$ ; Lens redshift
            strcompress(string(bcg_ra - lens5[0]/(60d*cos((!dpi/180d)*bcg_dec)),format=realformat),/remove_all)+'  ',$ ; x_ctr position
            strcompress(string(bcg_dec+ lens5[1]/60d,format=realformat),/remove_all)+'  ',$ ; y_ctr position
            strcompress(string(lens5[2],format=realformat),/remove_all)+'  ',$ ; Velocity dispersion of the NIE in 10^3 km/s
            strcompress(string(lens5[3],format=realformat),/remove_all)+'  ',$ ; constant kappa offset
            strcompress(string(lens5[4],format=realformat),/remove_all)+'  ',$ ; core radius (in arcmin)
            strcompress(string(lens5[5],format=realformat),/remove_all)+'  ',$ ; ellipticity (a-b)/(a+b)
            strcompress(string(lens5[6],format=realformat),/remove_all)        ; ellipse position angle
     

; Put in the regularization information
     printf,wunit,'# Regularization parameters for '
     printf,wunit,'#    kappa, gamma1, gamma2, alpha1, alpha2, flexion'
     printf,wunit,$
            strcompress(string(1d2,format=realformat),/remove_all)+'  ',$
            strcompress(string(0d,format=realformat),/remove_all)+'  ',$
            strcompress(string(0d,format=realformat),/remove_all)+'  ',$
            strcompress(string(0d,format=realformat),/remove_all)+'  ',$
            strcompress(string(0d,format=realformat),/remove_all)+'  ',$
            strcompress(string(0d,format=realformat),/remove_all)

     printf,wunit,'# Field parameters and BCG position'
     printf,wunit,$
            strcompress(string(field[0],format=realformat),/remove_all)+'  ',$ ; Field limits relative to the BCG
            strcompress(string(field[1],format=realformat),/remove_all)+'  ',$
            strcompress(string(field[0],format=realformat),/remove_all)+'  ',$
            strcompress(string(field[1],format=realformat),/remove_all)
     printf,wunit,$
            strcompress(string(bcg_ra,format=realformat),/remove_all)+'  ',$
            strcompress(string(bcg_dec,format=realformat),/remove_all)

     printf,wunit,'# Cluster member data'
     printf,wunit,'wcs'
     printf,wunit,strcompress(string(0L),/remove_all) ; Number of cluster members
; The rest of the CM data would have to go here...
     
     close,wunit
     free_lun,wunit


;; Put in a refinement file
     openw,wunit,name+tag+'_refine',/get_lun,width=18*10

     r2=1.2d
     r4=0.6d
     rs=0.15d

     printf,wunit,'# If the line starts with "0 0 radius refinement_level" '
     printf,wunit,'# we will refine arround each multiple image that is used'
     printf,wunit,'# The other regions are given by "RA Dec radius refinment_level" '
     printf,wunit,$
            strcompress(string(bcg_ra,format=realformat),/remove_all)+'  ',$
            strcompress(string(bcg_dec,format=realformat),/remove_all)+'  ',$
            strcompress(string(r2,format=realformat),/remove_all)+'  ',$
            strcompress(string(2,format=int_format),/remove_all)
     printf,wunit,$
            strcompress(string(bcg_ra,format=realformat),/remove_all)+'  ',$
            strcompress(string(bcg_dec,format=realformat),/remove_all)+'  ',$
            strcompress(string(r4,format=realformat),/remove_all)+'  ',$
            strcompress(string(4,format=int_format),/remove_all)
     printf,wunit,$
            strcompress(string(0d,format=realformat),/remove_all)+'  ',$
            strcompress(string(0d,format=realformat),/remove_all)+'  ',$
            strcompress(string(rs,format=realformat),/remove_all)+'  ',$
            strcompress(string(4,format=int_format),/remove_all)
     close,wunit
     free_lun,wunit


; Make the starting file for SWUnitedAMR
     openw,wunit,'run.'+name+tag,/get_lun,width=18*10
     printf,wunit,name+tag+'_lensmodel'
     printf,wunit,'0 '+'../'+file_basename(header)
     printf,wunit,name+tag+'_wl.cat'
     printf,wunit,'0 '+'../'+file_basename(header)
     printf,wunit,name+tag+'_sl.cat'
     printf,wunit,'0 '+'../'+file_basename(header)
     printf,wunit,name+tag+'_fl.cat'
     close,wunit
     free_lun,wunit


; Create the noise for the data
     
; WL noise is ellipticity
     e_src=sig_e*dcomplex(randomn(seed,nsrc,/double),randomn(seed,nsrc,/double))*nflag
     gt1=where(abs(e_src) gt 1d,ngt1)

; FL noise is intrinsic flexion
     intflex1 = sig_flex1*dcomplex(randomn(seed,nsrc,/double),randomn(seed,nsrc,/double))*nflag
     intflex3 = sig_flex3*dcomplex(randomn(seed,nsrc,/double),randomn(seed,nsrc,/double))*nflag

; SL noise is position errors
     dtheta=sig_theta*randomn(seed,nsrc,2,/double)*nflag

; Build the observed fields

; Add position noise
     theta_obs=cat[*,2:3] + dtheta ; in arcmin
     RADec_obs=cat[*,0:1]*0d       ; in degrees
     RADec_obs[*,0]=bcg_ra  - (theta_obs[*,0])/(60d*cos((!dpi/180d)*bcg_dec)) 
     RADec_obs[*,1]=bcg_dec + (theta_obs[*,1])/60d

; Convert source ellipiticity and shear into observed ellipticity
     gt1=where(abs(g) gt 1d,ngt1,complement=le1,ncomplement=nle1)
     e_img=e_src*0d
     if nle1 gt 0 then e_img[le1]=(e_src[le1] + g[le1])/(dcomplex(1d,0d) + conj(g[le1])*e_src[le1])
     if ngt1 gt 0 then e_img[gt1]=(dcomplex(1d,0d) + conj(e_src[gt1])*g[gt1])/conj(e_src[gt1] + g[gt1])
  
; Convert flexion and flexion noise into observed flexion
     psi1=flex1 + intflex1
     psi3=flex3 + intflex3

; Split up the images into SL, WL and FL measurements
     mu=1d/((1d - cat[*,13])^2 - total(cat[*,14:15]^2,2))

     multobj=where(cat[*,21] gt 1,n_mult,complement=singobj,ncomplement=n_sing)

     wlobj= where((cat[*,21] eq 1) and $
                  (abs(mu) gt mulim) and $
                  (abs(e_img) lt shearlim), n_wl)

     flobj= where((cat[*,21] eq 1) and $
                  (abs(mu) gt mulim) and $
                  (abs(psi1) lt flexmax1) and (abs(psi3) lt flexmax3) and $
                  (abs(psi1) gt flexmin1) and (abs(psi3) gt flexmin3), n_fl)

     slobj= where((cat[*,21] gt 1) and $
                  (abs(mu) gt mulim), n_sl)


; Calculate the chi-squared values for the true model vs the data
     eobs=e_img[wlobj]
     ebar=g[wlobj]
     big=where(abs(ebar) gt 1d,nbig)
     if nbig gt 0 then ebar[big]=1d/conj(ebar[big])
     chisq_wl=total(abs(eobs - ebar)^2/sig_e^2)
     chisq_f1=total(abs(intflex1[flobj]/sig_flex1)^2)
     chisq_f3=total(abs(intflex3[flobj]/sig_flex3)^2)
     
; Now the SL images, which is a bit more complicated
     b_sl=theta_obs[slobj,*] - my_lens(theta_obs[slobj,*],lens,z_weight=cat[slobj,5],/alpha_only)
     fct=1d/abs(mu[slobj])
     big=where(fct gt 1d3,nb)
     wee=where(fct lt 1d-4,nw)

     if nb gt 0 then fct[big]=1d3
     if nw gt 0 then fct[wee]=1d-4

     chisq_sl=total((b_sl - cat[slobj,8:9])^2/(sl_err*fct)^2)
     printf,notesunit,'Realization',realization
     printf,notesunit,'n_wl',n_wl
     printf,notesunit,'n_sl',n_sl
     printf,notesunit,'n_fl',n_fl
     printf,notesunit,'WL  ',chisq_wl
     printf,notesunit,'SL  ',chisq_sl
     printf,notesunit,'FL-1',chisq_f1
     printf,notesunit,'FL-3',chisq_f3
     
     print,'Realization',realization
     print,'Weak   ',n_wl
     print,'Strong ',n_sl
     print,'Flexion',n_fl
     print,''

     if (n_wl gt 0) then begin
; Save the WL data in the SWunited format
        wlout=dblarr(n_wl,6)
        wlout[*,0]=cat[wlobj,0]            ; RA
        wlout[*,1]=cat[wlobj,1]            ; Dec
        wlout[*,2]=real_part(e_img[wlobj]) ; e1
        wlout[*,3]=imaginary(e_img[wlobj]) ; e2
        wlout[*,4]=1d                      ; relative weight
        wlout[*,5]=cat[wlobj,4]            ; redshift


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
        sl_mu=mu[slobj]
        sl_theta=RADec_obs[slobj,*]

; Count the number of systems and the number of images in each system.
        nsys=n_elements(uniq(sl_cat[*,20]))
        nimg=sl_cat[uniq(sl_cat[*,20]),21]
        sys0=total(nimg,/cumulative)-nimg ; The indices in sl_cat of the first image in each system

        order=reverse(sort(sl_cat[sys0,21])) ; Sort them by decreasing # of imgs
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
     
        z_sys=sl_cat[sys0,4]
        
        for i=0,nsys-1 do begin
           ims=sys0[i]+lindgen(nimg[i]) ; Indices to the images in this family
           usesys=1

           printf,wunit,'# Image Family '+strcompress(string(sl_cat[sys0[i],20],format=int_format),/remove_all)
           printf,wunit,$
                  strcompress(string(i,format=int_format),/remove_all)+'  ',$        ; System index
                  strcompress(string(nimg[i],format=int_format),/remove_all)+'  ',$  ; Number of images
                  strcompress(string(usesys,format=int_format),/remove_all)+'  ',$   ; Whether to use
                  strcompress(string(z_sys[i],format=realformat),/remove_all)+'  ',$ ; System redshift
                  strcompress(string(0.01d,format=realformat),/remove_all)           ; Redshift error

           for j=0,nimg[i]-1 do begin
              m=ims[j]
              printf,wunit,$
                     strcompress(string(sl_theta[m,0],format=realformat),/remove_all)+'  ',$ ; Position
                     strcompress(string(sl_theta[m,1],format=realformat),/remove_all)+'  ',$
                     strcompress(string(sl_err,format=realformat),/remove_all)+'  ',$             ; Position errors
                     strcompress(string(sl_err,format=realformat),/remove_all)+'  ',$             ; (in decimal degrees)
                     strcompress(string(real_part(sl_e[m]),format=realformat),/remove_all)+'  ',$ ; Ellipticity
                     strcompress(string(imaginary(sl_e[m]),format=realformat),/remove_all)+'  ',$
                     strcompress(string(0.01d,format=realformat),/remove_all)+'  ',$ ; Ellipticity errors
                     strcompress(string(0.01d,format=realformat),/remove_all)+'  ',$ 
                     strcompress(string(sl_mu[m],format=realformat),/remove_all)+'  ',$ ; Flux - isn't used, so put in mu for ref.
                     strcompress(string(0.1d,format=realformat),/remove_all)            ; Flux error
           endfor
           
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
        flout[*,9]=cat[flobj,4]
     

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

  endfor


; Save a catalog of the observed quantities
  obscat=dblarr(nsrc,15)
  obscat[*,0:1]=RADec_obs
  obscat[*,2:3]=theta_obs
  obscat[*,4:5]=cat[*,4:5]
  obscat[*,6]= real_part(e_img)
  obscat[*,7]= imaginary(e_img)
  obscat[*,8]= real_part(psi1)
  obscat[*,9]= imaginary(psi1)
  obscat[*,10]=real_part(psi3)
  obscat[*,11]=imaginary(psi3)
  obscat[*,12:14]=cat[*,20:22]
  save_data,obscat,name+'_obs.dat',comment=' RA, Dec, theta1, theta2, redshift, Z(z), '+$
            'e1, e2, Psi11, Psi12, Psi31, Psi32, Image_family, N_images, i_image'

  close,notesunit
  free_lun,notesunit

  print,'Catalogs made!'
    
end
