function calc_moments, IMAGE, X_ORDER, Y_ORDER

; This function returns moments calculated about the center of the
; image.  There is no normalization, so from a "units" standpoint, the
; units are [LENGTH]^(X_ORDER+Y_ORDER).  If either X_ORDER or Y_ORDER
; is less than zero, it will be set to zero.  All moments are returned
; to double precision.  Furthermore, if the sum of X_ORDER and Y_ORDER
; is greater than unity (meaning second order or higher), the moments
; will be centroid-relative

; Find the dimensions and the center of the image
dims=double(size(image,/dimensions))
ctr=(dims-1d)/2d

; Make the coordinate arrays
x=(dindgen(dims) mod dims[0])-ctr[0]
y=double(lindgen(dims)/long(dims[0]))-ctr[1]

if x_order lt 0d then x_order=0d
if y_order lt 0d then y_order=0d

if x_order+y_order gt 1d then begin
   cts=calc_moments(image,0,0)
   ctr=[calc_moments(image,1,0),calc_moments(image,0,1)]/cts
endif else ctr=[0d,0d]

moment=total(image*((x-ctr[0])^x_order)*((y-ctr[1])^y_order),/double)

return,moment

end
