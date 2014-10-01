function invert_acs_dist, DIST_POS, CHIP1=CHIP1,CHIP2=CHIP2,$
                          F775W=F775W,F475W=F475W,NITER=NITER

if n_elements(dist_pos) ne 2 then begin
   print,'INVERT_ACS_DIST.PRO Error!'
   print,'This function only inverts one position at a time.'
   print,'Returning input!!'
   return,dist_pos
endif

; Inputs to the fitting
;fname='acs_dist_coord_deviate'
;fargs={dist_pos:dist_pos, chip1:keyword_set(chip1), $
;       chip2:keyword_set(chip2), f475w:keyword_set(f475w),$
;       f775w:keyword_set(f775w)}

; Set parameter limits
;parinf=replicate({value:0d, fixed:0, limited:[1,1], limits:[0d,0],$
;                  parname:''},2)

; Start at the center of the chip
;parinf[0].value=2048d
;parinf[1].value=1024d

; Keep it on the chip
;parinf[0].limits=[0d,4096d]
;parinf[1].limits=[0d,1024d]

;parinf[0].parname='X_DETECTOR'
;parinf[1].parname='Y_DETECTOR'

;detector_pos=mpfit(fname,functargs=fargs,bestnorm=chisq,dof=dof)


; We'll use a gradient method instead
;Start at the center of the chip.
det_pos_now=[2048d,1024d]

nmax=500
tolerance=1d

dcoo=0.1d
xhat=[1d,0d]
yhat=[0d,1d]

ssize=5d

for i=0,nmax-1 do begin

   dDISTANCE_dXdet=$
      (sqrt(total((dist_pos-$
                   acs_geom_distortion(det_pos_now+dcoo*xhat,$
                                       chip1=keyword_set(chip1),$
                                       chip2=keyword_set(chip2),$
                                       f775w=keyword_set(f775w),$
                                       f475w=keyword_set(f475w)))^2))-$
       sqrt(total((dist_pos-$
                   acs_geom_distortion(det_pos_now,$
                                       chip1=keyword_set(chip1),$
                                       chip2=keyword_set(chip2),$
                                       f775w=keyword_set(f775w),$
                                       f475w=keyword_set(f475w)))^2)))/dcoo
   
   dDISTANCE_dYdet=$
      (sqrt(total((dist_pos-$
                   acs_geom_distortion(det_pos_now+dcoo*yhat,$
                                       chip1=keyword_set(chip1),$
                                       chip2=keyword_set(chip2),$
                                       f775w=keyword_set(f775w),$
                                       f475w=keyword_set(f475w)))^2))-$
       sqrt(total((dist_pos-$
                   acs_geom_distortion(det_pos_now,$
                                       chip1=keyword_set(chip1),$
                                       chip2=keyword_set(chip2),$
                                       f775w=keyword_set(f775w),$
                                       f475w=keyword_set(f475w)))^2)))/dcoo

   gradient=[dDISTANCE_dXdet,dDISTANCE_dYdet]

   det_pos_new=det_pos_now - ssize*gradient

   DISTANCE_new=$
      sqrt(total((dist_pos-$
                  acs_geom_distortion(det_pos_new,$
                                      chip1=keyword_set(chip1),$
                                      chip2=keyword_set(chip2),$
                                      f775w=keyword_set(f775w),$
                                      f475w=keyword_set(f475w)))^2))
   if distance_new le tolerance then begin
      niter=i+1
      return,det_pos_new 
   endif else det_pos_now=det_pos_new

endfor

print,'ACS distortion inversion did not converge in the allotted iterations!'
return,det_pos_now

end
