FUNCTION BLOBMAKER, IMAGE, X0, MIN=MIN, NBLOB=NBLOB

; This function creates a list of the pixels in the 2D array IMAGE
; which are greater than MIN and are contiguous with X0.  If MIN is
; not specified, then it returns all positive pixels.  If IMAGE[X0] is
; not greater than MIN, then BLOBMAKER returns -1.

if not keyword_set(min) then min=0d

nblob=0
blob=-1

imdims=size(image,/dimensions)

if n_elements(x0) gt 1 then x=(lindgen(imdims))[x0[0],x0[1]] else x=x0

if image[x] lt min then return,blob

new=x
nnew=1
blob=new
while nnew gt 0 do begin

; Make a list of all pixels that adjoin a new pixel.
   test=-1
   for i=0,nnew-1 do begin
      newxy=array_indices(image,new[i])

      if newxy[0]+1 le imdims[0]-1 then $
         test=set_union(test,(lindgen(imdims))[newxy[0]+1,newxy[1]])
      if newxy[1]+1 le imdims[1]-1 then $
         test=set_union(test,(lindgen(imdims))[newxy[0],newxy[1]+1])
      if newxy[0]-1 ge 0 then $
         test=set_union(test,(lindgen(imdims))[newxy[0]-1,newxy[1]])
      if newxy[1]-1 ge 0 then $
         test=set_union(test,(lindgen(imdims))[newxy[0],newxy[1]-1])
   endfor

; Remove any pixels already in the blob
   test=set_difference(test,blob,count=ntest)

; If there are pixels to test, then test them.
   if ntest gt 0 then new=where(image[test] gt min,nnew) else begin
      new=-1
      nnew=0
   endelse

; The pixels that pass the test (if any) are added to the blob.
   if nnew gt 0 then new=test[new]
   blob=set_union(blob,new,count=nblob)

endwhile

RETURN,BLOB

END
