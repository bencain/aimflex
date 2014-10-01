function tinytim_psf, input_catalog, galaxies, psftag

; This will facilitate using TinyTim to create a set of PSF images.

xs=input_catalog[galaxies,1]
ys=input_catalog[galaxies,2]

chip1gals=where(ys gt (2052.886d +((2195d - 2052.886d)/(4152d - 105.327d))*$
                       (xs - 105.327d)),$
                nc1,complement=chip2gals,ncomplement=nc2)

ttpath='/usr/local/tinytim7.0/src/'

if nc1 gt 0 then begin

   chip1_img_xy=input_catalog[galaxies[chip1gals],1:2]
   chip1_det_xy=chip1_img_xy*0d

   for i=0,nc1-1 do begin
      forloop_status,i,nc1,wnum,$
                     label='Chip1: Inverting the ACS geometric distortion'
      chip1_det_xy[i,*]=invert_acs_dist(chip1_img_xy[i,*],$
                                        /chip1,/f775w)
   endfor
   forloop_status,0,0,wnum,/delete
;   delvarx,wnum
   save_data,chip1_det_xy,psftag+'_chip1.lis',/noheader

   window,wnum,xsize=800,ysize=500,/free,title='Correct tiny1 responses:'

   xyouts,0.05,0.8,'Camera:            15 (ACS WFC)',/normal,charsize=2
   xyouts,0.05,0.7,'Detector:          1',/normal,charsize=2
   xyouts,0.05,0.6,'Position:          @'+psftag+'_chip1.lis',/normal,$
          charsize=2
   xyouts,0.05,0.5,'Filter:            MONO',/normal,charsize=2
   xyouts,0.05,0.4,'Wavelength (nm):   775 (or the imaging filter)',/normal,$
          charsize=2
   xyouts,0.05,0.3,'PSF Size (arcsec): 1.5',/normal,charsize=2
   xyouts,0.05,0.2,'PSF file rootname: '+psftag+'_chip1img',/normal,$
          charsize=2
   
   spawn,ttpath+'tiny1 '+psftag+'_chip1.par'
   
   wdelete

   spawn,ttpath+'tiny2 '+psftag+'_chip1.par'
   for i=0,nc1-1 do spawn,ttpath+'tiny3 '+psftag+'_chip1.par pos='+$
                          strcompress(string(i),/remove_all)
endif

if nc2 gt 0 then begin

   chip2_img_xy=input_catalog[galaxies[chip2gals],1:2]
   chip2_det_xy=chip2_img_xy*0d

   for i=0,nc2-1 do begin
      forloop_status,i,nc2,wnum,$
                     label='Chip2: Inverting the ACS geometric distortion'
      chip2_det_xy[i,*]=invert_acs_dist(chip2_img_xy[i,*],$
                                        /chip2,/f775w)
   endfor
   forloop_status,0,0,wnum,/delete
   delvarx,wnum
   save_data,chip2_det_xy,psftag+'_chip2.lis',/noheader

   window,xsize=800,ysize=500,/free,title='Correct tiny1 responses:'

   xyouts,0.05,0.8,'Camera:            15 (ACS WFC)',/normal,charsize=2
   xyouts,0.05,0.7,'Detector:          2',/normal,charsize=2
   xyouts,0.05,0.6,'Position:          @'+psftag+'_chip2.lis',/normal,$
          charsize=2
   xyouts,0.05,0.5,'Filter:            MONO',/normal,charsize=2
   xyouts,0.05,0.4,'Wavelength (nm):   775 (the imaging filter)',/normal,$
          charsize=2
   xyouts,0.05,0.3,'PSF Size (arcsec): 1.5',/normal,charsize=2
   xyouts,0.05,0.2,'PSF file rootname: '+psftag+'_chip2img',/normal,$
          charsize=2

   spawn,ttpath+'tiny1 '+psftag+'_chip2.par'

   wdelete

   spawn,ttpath+'tiny2 '+psftag+'_chip2.par'
   for i=0,nc2-1 do spawn,ttpath+'tiny3 '+psftag+'_chip2.par pos='+$
                          strcompress(string(i),/remove_all)
endif

psf_files=strarr(n_elements(galaxies))

for i=0, nc1-1 do begin
   n_digits=floor(alog10(nc1)+1d,/L64)
   i_digits=floor(alog10(max([1L,i]))+1d,/L64)
;   print,i,n_digits,i_digits

   if n_digits eq 1 then fill='0' else begin
      fill=''
      for j=0,n_digits-i_digits-1 do fill+='0'
   endelse

   psf_files[chip1gals[i]]=psftag+'_chip1img'+fill+$
                           strcompress(string(i),/remove_all)+'.fits'
endfor

for i=0, nc2-1 do begin
   n_digits=floor(alog10(nc2)+1d,/L64)
   i_digits=floor(alog10(max([1L,i]))+1d,/L64)
;   print,i,n_digits,i_digits

   if n_digits eq 1 then fill='0' else begin
      fill=''
      for j=0,n_digits-i_digits-1 do fill+='0'
   endelse

   psf_files[chip2gals[i]]=psftag+'_chip2img'+fill+$
                           strcompress(string(i),/remove_all)+'.fits'
endfor

return,psf_files

end
