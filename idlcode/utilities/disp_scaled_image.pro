pro disp_scaled_image, img, wnum, delete=delete,title=title,$
                       label=label,winsize=winsize,log=log

if not keyword_set(title) then title='*** Image Display ***'
if not keyword_set(label) then label=''
if not keyword_set(winsize) then winsize=450

if ((n_elements(wnum) eq 0) or  (!d.window eq -1)) then begin
   window,xsize=winsize,ysize=winsize,/free,title=title
   wnum=!d.window
endif else begin
   wset,wnum
endelse


if keyword_set(delete) then wdelete,wnum else begin

   if keyword_set(log) then begin
      min=min(img)
      useimg=alog(img-min+1d)
      min=min(useimg)
      max=max(useimg)

      sizescale=double(winsize)/max(size(useimg,/dimensions))
      normscale=255d/max([(max-min),1d-6*abs(max)])
      
      tv,scale_image(normscale*(useimg-min),sizescale)

      xyouts,0.25,0.25,label,/normal,charsize=3.0
   endif else begin
      max=max(img)
      min=min(img)
   
      sizescale=double(winsize)/max(size(img,/dimensions))
      normscale=255d/max([(max-min),1d-6*abs(max)])
      
      tv,scale_image(normscale*(img-min),sizescale)

      xyouts,0.25,0.25,label,/normal,charsize=3.0
   endelse 
endelse


end
