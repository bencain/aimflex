function mk_window, size, radius=radius

; Make the size an integer
  size=floor(size,/l64)

; Set the size to just a 2-element array, plus some error control.
  if n_elements(size) ne 2 then size=[size[0],size[0]]

; We want an odd number for size 
  size=size+((size+1) mod 2L)

; Mark the midpoint of the array
  ctr=(size-1)/2
  
; Set the coordinates
  x=double(lindgen(size) mod size[0]-ctr[0])
  y=double(lindgen(size)/long(size[0])-ctr[1])

; Deal with the radius
  if not keyword_set(radius) then radius=double(min(ctr))-1d

; Make the window
  window=double((x^2+y^2) le radius^2)

  return, window
end

