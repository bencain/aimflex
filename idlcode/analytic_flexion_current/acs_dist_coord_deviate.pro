function acs_dist_coord_deviate, DET_POS, DIST_POS=DIST_POS,$
                        CHIP1=CHIP1,CHIP2=CHIP2,$
                        F775W=F775W,F475W=F475W

dev=dist_pos - acs_geom_distortion(det_pos,$
                                   chip1=keyword_set(chip1),$
                                   chip2=keyword_set(chip2),$
                                   f775w=keyword_set(f775w),$
                                   f475w=keyword_set(f475w))
;print,dev
return,dev/5d

end
