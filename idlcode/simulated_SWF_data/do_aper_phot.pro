pro do_aper_phot,fitsname,c1,c2,r,z_lens,fluxfile,wcs=wcs,clobber=clobber,tag=tag,buffer=buffer

  if n_elements(tag) gt 0 then tag=' '+strcompress(string(tag),/remove_all) else tag=''
  if n_elements(r) lt n_elements(c1) then r=replicate(r[0],n_elements(c1))
  if not keyword_set(buffer) then buffer=10d else buffer=buffer[0] ; background estimation annulus width

  realformat='(E15.8)'
  intformat='(I6)'

  image=mrdfits(fitsname,/silent)
  hdr=headfits(fitsname)

  dims=size(image,/dimension)
  imx=(dindgen(dims) mod dims[0])+1d ;Note that these are FITS image pixel locations
  imy=double(lindgen(dims)/long(dims[0]))+1d

; Set NANs to zero
  nans=where(finite(image) eq 0,n_nans)
  if n_nans gt 0 then begin
     print,string(n_nans,format=intformat)+' bad pixel values detected - setting them to zero...'
     image[nans]=0d
  endif

  if not keyword_set(clobber) then begin
     outexists=file_test(fluxfile)
     if outexists then begin
;        print,'updating existing file'
        openu,out,fluxfile,/get_lun,/append,width=200 
     endif else begin
;        print,'opening new file'
        openw,out,fluxfile,/get_lun,width=200
     endelse
  endif else begin
;     print,'clobbering existing file'
     outexists=0
     openw,out,fluxfile,/get_lun,width=200
  endelse

;;;;;;;;;;;
  extast,hdr,astr,np
  if np ne 2 then exit
  cd = astr.cd
  cdelt = astr.cdelt
  crval = astr.crval
  if cdelt[0] NE 1.0 then begin
     cd[0,0] *= cdelt[0] 
     cd[0,1] *= cdelt[0]
     cd[1,1] *= cdelt[1]
     cd[1,0] *= cdelt[1]
  endif

  cdinv = invert(cd,/double)
  crpix = astr.crpix
;;;;;;;;;;;

; Create WCS pixel positions
  imdx=imx-crpix[0]
  imdy=imy-crpix[1]
  
  ima = cd[0,0]*imdx + cd[0,1]*imdy ;Can't use matrix multiplication, imdx and imdy might be vectors
  imd = cd[1,0]*imdx + cd[1,1]*imdy
 
  ima += crval[0]
  imd += crval[1]



  if keyword_set(wcs) then begin
     ; Convert from WCS coordinates to pixel coordinates
     print,'Starting from WCS Coords'
     a=c1
     d=c2

     ad2xy,a,d,astr,x,y

     pscale=sqrt(abs(determ(astr.cd)))*3600d ; plate scale in arcsec/pixel
     rwcs=r
     rpix=rwcs/pscale

  endif else begin
     print,'Starting from Image Coords'
     x=c1
     y=c2

     xy2ad,x,y,astr,a,d
     
     pscale=sqrt(abs(determ(astr.cd)))*3600d ; plate scale in arcsec/pixel
     rpix=r
     rwcs=rpix*pscale

  endelse


  flux=0d*a
  nok=0L*long(a)
  bkg_per_pix=flux


  for i=0,n_elements(a)-1 do begin

     okpix=where(sqrt((imx-x[i])^2+(imy-y[i])^2) le rpix[i],n)     
     nok[i]=n
;     if nok[i] gt 0L then flux[i]=total(image[okpix]) else flux[i]=0d ; This is in units of pixel^2
;;;;;;;;;;
;     !p.multi=[0,2,1]
     if nok[i] gt 0L then begin
        bkgpix=where((sqrt((imx-x[i])^2+(imy-y[i])^2) gt rpix[i]) and (sqrt((imx-x[i])^2+(imy-y[i])^2) le rpix[i]+buffer),n)
        
;        bkgdata=dblarr(3,n)
;        bkgdata[0,*]=imx[bkgpix]
;        bkgdata[1,*]=imy[bkgpix]
;        bkgdata[2,*]=image[bkgpix]

;        bkgfit=sfit(bkgdata,1,/irregular,kx=coeff,/max_degree)

;        bkginterp=coeff[0]+$
 ;                 coeff[1]*imy[okpix]+$
  ;                coeff[2]*imx[okpix]

;        bkginterp=coeff[0]+$
 ;                 coeff[1]*imy[okpix]+$
  ;                coeff[2]*imy[okpix]^2+$
   ;               coeff[3]*imx[okpix]+$
    ;              coeff[4]*imx[okpix]*imy[okpix]+$
     ;             coeff[5]*imx[okpix]^2
      ;  bkginterp=mean(image[bkgpix])

;        bkg_sub_img=image[okpix] - bkginterp

;        bkg_sub_img=image[okpix] - mean(image[bkgpix])


;        plot,imx[okpix],image[okpix],/psym,charsize=2
;        oplot,imx[okpix],bkginterp,psym=4,color=fsc_color('green')

;        plot,imy[okpix],image[okpix],/psym,charsize=2
;        oplot,imy[okpix],bkginterp,psym=4,color=fsc_color('green')
 
;        cgcontour,bkg_sub_img,imx[okpix],imy[okpix],/irregular,/isotropic
        
;        wait,1.5


        flux[i]=total(image[okpix])
        bkg_per_pix[i]=min(image[bkgpix])
        
        

     endif else flux[i]=0d      ; This is in units of pixel^2

;;;;;;;;;;
  endfor


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  c=2.99792458d5                ; km/s
  G=4.302d-9                    ; (km/s)^2(Mpc/Msun)
  H0=70d                        ; km/s/Mpc
  DH=c/H0                       ; Mpc

  z_src=1d4
  DL=redshift_to_angdd(z_lens)*DH          ; Mpc
  DS=redshift_to_angdd(z_src)*DH           ; Mpc
  DLS=DS - ((1d + z_lens)/(1d + z_src))*DL ; Mpc

  Sig_c=(c^2/(4d*!dpi*G))*(DS/(DL*DLS)) ; Msun/Mpc^2
  fct=Sig_c*DL^2                        ; Now in Msun/Sr
  fct*=(!dpi/(180d*3600d))^2            ; Now in Msun/arcsec^2
  fct*=pscale^2                         ; Now in Msun/pixel^2



  if not outexists then printf,out,'# MFile region x[pix] y[pix] RA[deg] Dec[deg] r[pix] r["] Mass[Msun] NAperPix BkgPerPixel[Msun] Scale["/pix] DetFile'
  for i=0,n_elements(a)-1 do printf,out,$
                                    file_basename(fitsname)+' '+$
                                    string(i,format=intformat)+' '+$
                                    string(x[i],format=realformat)+' '+$
                                    string(y[i],format=realformat)+' '+$
                                    string(a[i],format=realformat)+' '+$
                                    string(d[i],format=realformat)+' '+$
                                    string(rpix[i],format=realformat)+' '+$
                                    string(rwcs[i],format=realformat)+' '+$
                                    string(flux[i]*fct,format=realformat)+' '+$
                                    string(nok[i],format=intformat)+' '+$
                                    string(bkg_per_pix[i]*fct,format=realformat)+ $
                                    string(pscale,format=realformat)+ $
                                    tag

  close,out
  free_lun,out
end
