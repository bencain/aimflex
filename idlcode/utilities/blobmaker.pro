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

if x ge n_elements(image) then return,blob
if x lt 0 then return,blob
if image[x] lt min then return,blob

labels=label_region(image gt min,/ulong)

blob=where(labels eq labels[x],nblob)

RETURN,BLOB

END
