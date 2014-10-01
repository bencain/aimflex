function clean_img, img, fitcat, fomcat, img_err, psf_fwhm=psf_fwhm,$
                    sersic=sersic

if keyword_set(sersic) then npar=13 else npar=12
if not keyword_set(psf_fwhm) then psf_fwhm=1.8d

; Cleaning a data image

; I assume that the catalogs input here are red-sequence galaxy images
good=where(fomcat[*,13] eq 0,ngood)
rad=fitcat[*,6]
newimg=img
psf=mk_gaussian_psf(psf_fwhm)

; Subtract out each of the red sequence galaxies
for i=0,ngood-1 do begin

   forloop_status,i,ngood,wnum
   
   stamprad=10d*rad[good[i]]
   stamprad=max([stamprad,sqrt(fomcat[*,2]/!dpi)])
   stamprad=min([stamprad,750d])

   win=mk_window(stamprad)
   fps=fitcat[good[i],3:npar+2]

   errors=2d*img_err*randomn(seed,size(win,/dimensions))*win

   neg_img=-aim_mk_image(fps,0d,win,psf,sersic=keyword_set(sersic))+errors

   newimg=insert_stamp(newimg,neg_img,fitcat[good[i],1:2],/sum)

endfor

forloop_status,0,0,wnum,/delete

return,newimg

end
