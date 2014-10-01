function insert_stamp, image, stamp, center, sum=sum

;This function inserts a stamp into a larger image.  The stamp is
;assumed to be have an odd number of pixels on each side and the
;existing image in the location of the stamp is replaced unless the
;SUM keyword is set, in which case the stamp is added to the image.
;The center of the stamp is aligned with the CENTER parameter.  The
;new image is returned.

newimg=image

center=long(center)

image_size=size(image,/dimensions)
stamp_size=size(stamp,/dimensions)
stamp_ctr=(stamp_size-1)/2L

; Find the corners of the stamp as projected onto the image
xmin=center[0]-stamp_ctr[0]
xmax=center[0]+stamp_ctr[0]

ymin=center[1]-stamp_ctr[1]
ymax=center[1]+stamp_ctr[1]

;print,xmin,xmax,ymin,ymax

; Find the overlaps on the sides
overlap_xup=max([xmax-image_size[0]+1,0])
overlap_xdn=max([0-xmin,0])

overlap_yup=max([ymax-image_size[1]+1,0])
overlap_ydn=max([0-ymin,0])

;print,overlap_xup,overlap_xdn,overlap_yup,overlap_ydn

; Insert the stamp
if keyword_set(sum) then begin
   newimg[(xmin+overlap_xdn):(xmax-overlap_xup),$
          (ymin+overlap_ydn):(ymax-overlap_yup)]+=$
      stamp[(overlap_xdn):(stamp_size[0]-overlap_xup-1),$
            (overlap_ydn):(stamp_size[1]-overlap_yup-1)]
endif else begin
;   print,(xmin+overlap_xdn),(xmax-overlap_xup),(ymin+overlap_ydn),(ymax-overlap_yup)
   newimg[(xmin+overlap_xdn):(xmax-overlap_xup),$
            (ymin+overlap_ydn):(ymax-overlap_yup)]=$
      stamp[(overlap_xdn):(stamp_size[0]-overlap_xup-1),$
            (overlap_ydn):(stamp_size[1]-overlap_yup-1)]
endelse

return,newimg

end
