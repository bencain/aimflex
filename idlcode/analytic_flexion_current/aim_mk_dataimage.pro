function aim_mk_dataimage, mpars, bkg, window, psf, $
                           nopsf=nopsf, $
                           noise=noise, gaussian_noise=gaussian_noise, $
                           seed=seed, $
                           floor=floor, $
                           sersic=sersic, pseudogauss=pseudogauss,$
                           moffat=moffat,$
                           parent_scale=parent_scale,$
                           centroid_tol=centroid_tol,$
                           replace_center=replace_center,$
                           max_iter=max_iter,$
                           parent_image=parent_image

if not keyword_set(parent_scale) then parent_scale=2d
if not keyword_set(centroid_tol) then centroid_tol=1.5d
if not keyword_set(max_iter) then max_iter=10

; This function uses AIM_MK_IMAGE to create a data image.  The main
; difference is that the returned data image has the centroid within
; CENTROID_TOL of the center of the window.

; Make the large window function for the parent image
wdims=size(window,/dimensions)
parent_window=mk_window(double(max(wdims))*parent_scale)

; Make the large parent image
parent_image=aim_mk_image(mpars,bkg,parent_window,psf,$
                         nopsf=keyword_set(nopsf),$
                         noise=keyword_set(noise),$
                         gaussian_noise=gaussian_noise,$
                         seed=seed,$
                         floor=keyword_set(floor),$
                         sersic=keyword_set(sersic),$
                         pseudogauss=keyword_set(pseudogauss),$
                         moffat=keyword_set(moffat))


; Initialize the centroid (we will keep this WRT the parent frame).
m10=calc_moments(parent_image-bkg*parent_window,1,0)
m01=calc_moments(parent_image-bkg*parent_window,0,1)
m00=calc_moments(parent_image-bkg*parent_window,0,0)

midpoint=(size(parent_window,/dimensions)-1d)/2d
centroid=[m10,m01]/m00

stamp=dblarr(size(window,/dimensions))

for i=0,max_iter-1 do begin
   stamp=cut_stamp(parent_image,centroid+midpoint,window)
   
   delta_centroid=[calc_moments(stamp-bkg*window,1,0),$
                   calc_moments(stamp-bkg*window,0,1)]/$
                  calc_moments(stamp-bkg*window,0,0)

   if sqrt(total(delta_centroid^2)) lt centroid_tol then break

   if sqrt(total((centroid+delta_centroid)^2)) gt min(midpoint) then break
   centroid+=delta_centroid
   
endfor

; We want to record the offset of the centroid.  We'll return
; -centroid

if keyword_set(replace_center) then mpars[1:2]=-centroid


return,stamp

end
