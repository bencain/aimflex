pro ubercat, name, density, fieldsize, header=header, noise=noise, bcg_ra=bcg_ra,bcg_dec=bcg_dec,$
             lens0=lens0, lens1=lens1, lens2=lens2, lens3=lens3, lens4=lens4, lens5=lens5, z_lens=z_lens, $
             grid=grid

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

; How many
  nsrc=ceil(density*fieldsize^2,/l64)

  multiplicity=10d              ; How many images per source to account for
  nim_max=nsrc*multiplicity

  cat=dblarr(nim_max,24)

; Draw the observed magnitude, the redshift, and the observed size
  mzr=get_mzr(nsrc,seed=seed,fixedm=24d,fixedr=0.3d)

  default=[0d,0d,$              ; Center
           0d,0d,1d,$           ; sigma, k0, theta_c
           0d,0d,0.94369d]       ; eps, phi, DLInf/DInf

  if not keyword_set(lens0) then lens0=default
  if not keyword_set(lens1) then lens1=default
  if not keyword_set(lens2) then lens2=default
  if not keyword_set(lens3) then lens3=default
  if not keyword_set(lens4) then lens4=default
  if not keyword_set(lens5) then lens5=default

  lens=[lens0,lens1,lens2,lens3,lens4,lens5]
  if not keyword_set(z_lens) then z_lens=0.2d

; Fitting tolerance
  err=1d-12

  count=0                       ; Number of images output
  used=0                        ; Number of sources used
  family=0


; If the option is set, find the betas that give you a gridded set of
; thetas with the same density.  We'll be adding in some
; multiple images too, but that's ok.
  if keyword_set(grid) then begin
     ng=long(ceil(sqrt(nsrc)))
     ts=dblarr(ng^2,2)
     ts[*,0]=(dindgen(ng^2) mod ng)/double(ng-1) - 0.5d
     ts[*,1]=double(lindgen(ng^2)/ng)/double(ng-1) - 0.5d
     ts*=fieldsize

; Now randomly permute the order of the sources
     mixed=(sort(randomu(seed,ng^2)))[0:nsrc-1]
     ts=ts[mixed,*]

     zs=dblarr(nsrc)
     for i=0,nsrc-1 do zs[i]=redshift_to_weight(mzr[i,1],z_lens)

     as=my_lens(ts,lens,z_weight=zs,/alpha_only)
     bs=ts - as

  endif


  while ((count lt nim_max) and (used lt nsrc)) do begin

     if ((used) mod 25) eq 0 then $
        print,"Finished "+strcompress(string(used,format=int_format),/remove_all)+$
              " of a possible "+strcompress(string(nsrc,format=int_format),/remove_all)+$
              " source positions..."

     if keyword_set(grid) then begin
        b=bs[used,*]
        z=zs[used]
     endif else begin
        b=fieldsize*(randomu(seed,2,/double)-0.5d)
        Z=redshift_to_weight(mzr[used,1],z_lens)
     endelse

     t=find_mult_images_triangles(b,err,z,lens,span=2d*fieldsize,ngrid=20L,nimages=nimages)

     start=count
     count+=nimages

     if ((nimages gt 0) and (count lt nim_max)) then begin

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

        lens_fields=my_lens(t,lens,z_weight=Z)

        mu=abs(1d/((1d - lens_fields[*,3])^2 - lens_fields[*,4]^2 - lens_fields[*,5]^2))

        cat[start:count-1, 0]=bcg_ra - t[*,0]/(60d*cos((!dpi/180d)*bcg_dec))      ; RA - observed
        cat[start:count-1, 1]=bcg_dec + t[*,1]/60d                                ; Dec -observed
        cat[start:count-1, 2]=t[*,0]                                              ; Theta1
        cat[start:count-1, 3]=t[*,1]                                              ; Theta2
        cat[start:count-1, 4]=mzr[used,0]+2.5d*alog10(mu)                         ; magnitude
        cat[start:count-1, 5]=mzr[used,2]*sqrt(mu)                                ; size, already in arcsec
        cat[start:count-1, 6]=mzr[used,1]                                         ; redshift
        cat[start:count-1, 7]=bcg_ra - b[0]/(60d*cos((!dpi/180d)*bcg_dec))        ; RA - unlensed
        cat[start:count-1, 8]=bcg_dec+ b[1]/60d                                   ; Dec -unlensed
        cat[start:count-1, 9]=b[0]                                                ; beta1
        cat[start:count-1,10]=b[1]                                                ; beta2
        cat[start:count-1,11]=lens_fields[*,0]                                    ; psi
        cat[start:count-1,12]=lens_fields[*,1]                                    ; alpha1
        cat[start:count-1,13]=lens_fields[*,2]                                    ; alpha2
        cat[start:count-1,14]=lens_fields[*,3]                                    ; kappa
        cat[start:count-1,15]=lens_fields[*,4]                                    ; gamma1
        cat[start:count-1,16]=lens_fields[*,5]                                    ; gamma2
        cat[start:count-1,17]=lens_fields[*,6]/60d                                ; fflex1   Note that the flexion is calculated
        cat[start:count-1,18]=lens_fields[*,7]/60d                                ; fflex2   in units of 1/arcmin, but we change
        cat[start:count-1,19]=lens_fields[*,8]/60d                                ; gflex1   that to 1/arcsec here.
        cat[start:count-1,20]=lens_fields[*,9]/60d                                ; gflex2
        cat[start:count-1,21]=family                                              ; Image family
        cat[start:count-1,22]=nimages                                             ; Number of images
        cat[start:count-1,23]=dindgen(nimages)                                    ; Image index

        family++
        used++
     endif

  endwhile

; Kick out foreground objects
  field=fieldsize*[-0.5d,0.5d]

  obj=where((cat[*,22] gt 0) and $
            (cat[*,6] gt z_lens), nobj)
  outcat=cat[obj,*]

; Center the field on the BCG (WCS coords in decimal degrees)
; Keep the left/right orientation correct w.r.t. RA, and account for
; the declination.

  save_data,outcat,name+'.dat',comment=' RA, Dec, theta1, theta2, magnitude, size, redshift, '+$
            'beta_RA, beta_Dec, beta1, beta2, psi, alpha1, alpha2, kappa, gamma1, gamma2, F1, F2, G1, G2, '+$
            'image_family, N_images, i_image'

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

; Plot the results of the deflections
  plot,outcat[*,2],outcat[*,3],psym=4,/isotropic
  oplot,outcat[*,9],outcat[*,10],psym=2,symsize=2,color=fsc_color('red')
  for i=0,nobj-1 do oplot,[outcat[i,2],outcat[i,9]],[outcat[i,3],outcat[i,10]],color=fsc_color('green')
  

;;;;;;;;;;;;

; Set up the other files...
; We need to have a lensmodel file...

  openw,wunit,name+'_lensmodel',/get_lun,width=18*10
  
  printf,wunit,'# Lens model parameter file, with BCG and other stuff also'
  printf,wunit,'wcs'
  
  printf,wunit,strcompress(string(1,format=int_format),/remove_all) ; 1 lens

;  DLI_over_DI=1d - ((1d + z_lens)*redshift_to_angdd(z_lens))/((1d + 1d4)*redshift_to_angdd(1d4))
  printf,wunit,$
         strcompress(string(2,format=int_format),/remove_all)+'  ',$        ; NIE lens
         strcompress(string(z_lens,format=realformat),/remove_all)+'  ',$   ; Lens redshift
         strcompress(string(bcg_ra,format=realformat),/remove_all)+'  ',$   ; x_ctr position
         strcompress(string(bcg_dec,format=realformat),/remove_all)+'  ',$  ; y_ctr position
         strcompress(string(lens0[2],format=realformat),/remove_all)+'  ',$ ; Velocity dispersion of the NIE in 10^3 km/s
         strcompress(string(lens0[3],format=realformat),/remove_all)+'  ',$ ; constant kappa offset
; We'll pretend that we're ignorant to any ellipticity or substructure.
         strcompress(string(lens0[4],format=realformat),/remove_all)+'  ',$    ; core radius (in arcmin)
         strcompress(string(lens0[5]*0d,format=realformat),/remove_all)+'  ',$ ; ellipticity (a-b)/(a+b)
         strcompress(string(lens0[6]*0d,format=realformat),/remove_all)        ; ellipse position angle

  printf,wunit,'# True model (main halo): '+$
         strcompress(string(2,format=int_format),/remove_all)+'  ',$        ; NIE lens
         strcompress(string(z_lens,format=realformat),/remove_all)+'  ',$   ; Lens redshift
         strcompress(string(bcg_ra,format=realformat),/remove_all)+'  ',$   ; x_ctr position
         strcompress(string(bcg_dec,format=realformat),/remove_all)+'  ',$  ; y_ctr position
         strcompress(string(lens0[2],format=realformat),/remove_all)+'  ',$ ; Velocity dispersion of the NIE in 10^3 km/s
         strcompress(string(lens0[3],format=realformat),/remove_all)+'  ',$ ; constant kappa offset
         strcompress(string(lens0[4],format=realformat),/remove_all)+'  ',$ ; core radius (in arcmin)
         strcompress(string(lens0[5],format=realformat),/remove_all)+'  ',$ ; ellipticity (a-b)/(a+b)
         strcompress(string(lens0[6],format=realformat),/remove_all)        ; ellipse position angle


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
  openw,wunit,name+'_refine',/get_lun,width=18*10

  printf,wunit,'# If the line starts with "0 0 radius refinement_level" '
  printf,wunit,'# we will refine arround each multiple image that is used'
  printf,wunit,'# The other regions are given by "RA Dec radius refinment_level" '
  printf,wunit,$
         strcompress(string(bcg_ra,format=realformat),/remove_all)+'  ',$
         strcompress(string(bcg_dec,format=realformat),/remove_all)+'  ',$
         strcompress(string(1.5d,format=realformat),/remove_all)+'  ',$
         strcompress(string(2,format=int_format),/remove_all)
  printf,wunit,$
         strcompress(string(bcg_ra,format=realformat),/remove_all)+'  ',$
         strcompress(string(bcg_dec,format=realformat),/remove_all)+'  ',$
         strcompress(string(2d/3d,format=realformat),/remove_all)+'  ',$
         strcompress(string(4,format=int_format),/remove_all)
  printf,wunit,$
         strcompress(string(0d,format=realformat),/remove_all)+'  ',$
         strcompress(string(0d,format=realformat),/remove_all)+'  ',$
         strcompress(string(0.1d,format=realformat),/remove_all)+'  ',$
         strcompress(string(4,format=int_format),/remove_all)
  close,wunit
  free_lun,wunit

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
  printf,notesunit,'nobj = ',nobj
  close,notesunit
  free_lun,notesunit
  


; Make the S, W and Flexion lensing catalogs from the full catalog
  split_cats,name,seed=seed,noise=keyword_set(noise), bcg_ra=bcg_ra, bcg_dec=bcg_dec

; Make the starting file for SWUnitedAMR
  openw,wunit,'run.'+name,/get_lun,width=18*10
  printf,wunit,name+'_lensmodel'
  printf,wunit,'0 '+'../'+file_basename(header)
  printf,wunit,name+'_wl.cat'
  printf,wunit,'0 '+'../'+file_basename(header)
  printf,wunit,name+'_sl.cat'
  printf,wunit,'0 '+'../'+file_basename(header)
  printf,wunit,name+'_fl.cat'
  close,wunit
  free_lun,wunit
  



  print,'Catalogs made!'




end
