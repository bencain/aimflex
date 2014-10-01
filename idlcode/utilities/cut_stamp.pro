function cut_stamp, parent_image, center, window, leave_zeroes=leave_zeroes

; Make sure that the center is a valid pixel.  x and y refer to the
; parent image. u and v will refer to the stamp
  x0=floor(center[0])
  y0=floor(center[1])

; Charaterize the stamp and the limits of the stamp on the parent
; image
  stampdims=size(window,/dimensions)
  parentdims=size(parent_image,/dimensions)

; Set the center point in the stamp.  
  u0=(stampdims[0]-1)/2
  v0=(stampdims[1]-1)/2

; Set the boundary of the overlap in the parent image coordinates.
  xlim=[max([0L,x0-u0]),min([parentdims[0]-1L,x0+u0])]
  ylim=[max([0L,y0-v0]),min([parentdims[1]-1L,y0+v0])]

; Set the boundary of the overlap in the stamp coordinates.
  ulim=[max([u0-x0,0L]),min([u0+parentdims[0]-x0-1L,2L*u0])]
  vlim=[max([v0-y0,0L]),min([v0+parentdims[1]-y0-1L,2L*v0])]

; Make the stamp and put in the parent values
  stamp=dblarr(stampdims)
  stamp[ulim[0]:ulim[1],$
        vlim[0]:vlim[1]]=parent_image[xlim[0]:xlim[1],$
                                      ylim[0]:ylim[1]]

; Window the stamp
  stamp*=window

; Cut out any overlapped portion of the window (e.g., if the stamp
; overhangs the edge of the parent image).
  if not keyword_set(leave_zeroes) then begin
     zeroes=where(stamp eq 0d,nzeroes)
     if nzeroes gt 0 then window[zeroes]=0d
  endif

  return, stamp

end
