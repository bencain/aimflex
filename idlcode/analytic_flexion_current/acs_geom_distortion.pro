function COORD_WARP,XY_POS,$
                    CHIP1=CHIP1,CHIP2=CHIP2,$
                    F775W=F775W,F475W=F475W

; Make sure that only one chip and one filter each are set.
if (keyword_set(chip1) and keyword_set(chip2)) or $
   ((not keyword_set(chip1)) and (not keyword_set(chip2))) then begin
   print,'** ACS_GEOM_DISTORTION.PRO ERROR**'
   print, 'Must set one (and only one) of CHIP1 and CHIP2 keywords.'
   return,xy_pos
endif

if (keyword_set(f775w) and keyword_set(f475w)) or $
   ((not keyword_set(f775w)) and (not keyword_set(f475w))) then begin
   print,'** ACS_GEOM_DISTORTION.PRO ERROR**'
   print, 'Must set one (and only one) of F475W and F775W keywords.'
   return,xy_pos
endif

; Now we'll be able to distort pixels from the detector frame
; to the sky frame.
cx=dblarr(5,5)
cy=dblarr(5,5)

; Reference centerpoint for each chip
xref=2048d
yref=1024d

if keyword_set(f475w) then begin

   if keyword_set(chip1) then begin

      cx[1,0]=2.3549034d-3
      cx[1,1]=4.9228288d-2
      cx[2,0]=8.7479350d-8
      cx[2,1]=-3.5462247d-7
      cx[2,2]=4.1167115d-7
      cx[3,0]=2.4919962d-12
      cx[3,1]=-2.6430174d-11
      cx[3,2]=-6.1309985d-12
      cx[3,3]=-2.3011334d-11
      cx[4,0]=6.8535452d-15
      cx[4,1]=-1.0477702d-15
      cx[4,2]=4.8720065d-15
      cx[4,3]=1.5629619d-15
      cx[4,4]=9.3663040d-16
      
      cy[1,0]=4.8579425d-2
      cy[1,1]=2.0384898d-3
      cy[2,0]=-4.7027800d-7
      cy[2,1]=2.9554738d-7
      cy[2,2]=-1.2650904d-7
      cy[3,0]=-1.8907148d-11
      cy[3,1]=-2.2915760d-12
      cy[3,2]=-2.2080882d-11
      cy[3,3]=3.5603606d-12
      cy[4,0]=-8.4046268d-15
      cy[4,1]=2.2479310d-15
      cy[4,2]=-4.8650557d-15
      cy[4,3]=-1.6374240d-15
      cy[4,4]=-8.1116375d-16
      
   endif else if keyword_set(chip2) then begin

      cx[1,0]=2.0538305E-03
      cx[1,1]=4.9831815E-02
      cx[2,0]=9.3801589E-08
      cx[2,1]=-2.4732915E-07
      cx[2,2]=4.2238801E-07
      cx[3,0]=-4.5194209E-13
      cx[3,1]=-2.6516291E-11
      cx[3,2]=-3.2822107E-12
      cx[3,3]=-2.3123720E-11
      cx[4,0]=-4.1951103E-17
      cx[4,1]=-7.8108164E-16
      cx[4,2]=1.7519748E-15
      cx[4,3]=5.0924065E-16
      cx[4,4]=1.1268242E-15
      
      cy[1,0]=5.0270893d-2
      cy[1,1]=1.3970418d-3
      cy[2,0]=-3.6028140d-7
      cy[2,1]=3.0452276d-7
      cy[2,2]=-7.6292842d-8
      cy[3,0]=-2.1354708d-11
      cy[3,1]=-4.5596170d-12
      cy[3,2]=-2.5823590d-11
      cy[3,3]=3.6116156d-12
      cy[4,0]=-1.0656479d-15
      cy[4,1]=-8.1105427d-16
      cy[4,2]=-1.2633409d-15
      cy[4,3]=-1.1973823d-16
      cy[4,4]=-7.5669360d-16

   endif
endif else if keyword_set(f775w) then begin
   if keyword_set(chip1) then begin

      cx[1,0]=2.3555558d-3
      cx[1,1]=4.9226720d-2
      cx[2,0]=8.7665057d-8
      cx[2,1]=-3.5492980d-7
      cx[2,2]=4.1204797d-7
      cx[3,0]=2.2064750d-12
      cx[3,1]=-2.6304434d-11
      cx[3,2]=-6.2123557d-12
      cx[3,3]=-2.2973637d-11
      cx[4,0]=6.6497368d-15
      cx[4,1]=-8.1341739d-16
      cx[4,2]=4.6849308d-15
      cx[4,3]=1.5693484d-15
      cx[4,4]=9.3017135d-16
      
      cy[1,0]=4.8577465d-2
      cy[1,1]=2.0394567d-3
      cy[2,0]=-4.7147256d-7
      cy[2,1]=2.9669084d-7
      cy[2,2]=-1.2672213d-7
      cy[3,0]=-1.9791956d-11
      cy[3,1]=-2.6351202d-12
      cy[3,2]=-2.1893152d-11
      cy[3,3]=3.3894151d-12
      cy[4,0]=-7.3565718d-15
      cy[4,1]=1.1993947d-15
      cy[4,2]=-5.0465719d-15
      cy[4,3]=-1.8030409d-15
      cy[4,4]=-8.3221601d-16

   endif else if keyword_set(chip2) then begin

      cx[1,0]=2.0536850d-3
      cx[1,1]=4.9831357d-2
      cx[2,0]=9.3125600d-8
      cx[2,1]=-2.4647176d-7
      cx[2,2]=4.2216712d-7
      cx[3,0]=-2.7969348d-13
      cx[3,1]=-2.6402314d-11
      cx[3,2]=-3.0637602d-12
      cx[3,3]=-2.3506011d-11
      cx[4,0]=5.5156753d-16
      cx[4,1]=-1.1282359d-15
      cx[4,2]=1.7506486d-15
      cx[4,3]=3.6400651d-16
      cx[4,4]=1.2064690d-15
      
      cy[1,0]=5.0269973d-2
      cy[1,1]=1.3971397d-3
      cy[2,0]=-3.6013776d-7
      cy[2,1]=3.0353459d-7
      cy[2,2]=-7.6010458d-8
      cy[3,0]=-2.1652093d-11
      cy[3,1]=-3.9016108d-12
      cy[3,2]=-2.6063049d-11
      cy[3,3]=3.6441852d-12
      cy[4,0]=-1.7951520d-16
      cy[4,1]=3.4434957d-16
      cy[4,2]=-1.3970546d-15
      cy[4,3]=-1.4914283d-16
      cy[4,4]=-8.1411815d-16

   endif
endif

x=xy_pos[*,0]
y=xy_pos[*,1]

xcorr=x*0d
ycorr=y*0d

for i=0,4 do for j=0,i do xcorr+=cx[i,j]*(x-xref)^j*(y-yref)^(i-j)
for i=0,4 do for j=0,i do ycorr+=cy[i,j]*(x-xref)^j*(y-yref)^(i-j)

; Return the warped coordinates
out_xy_pos=xy_pos
out_xy_pos[*,0]=xcorr
out_xy_pos[*,1]=ycorr

;print,'COORD_WARP',size(out_xy_pos,/dimensions)

return,out_xy_pos

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function coord_twochip,XY_POS,$
                     CHIP1=CHIP1,CHIP2=CHIP2,$
                     F775W=F775W,F475W=F475W

; Now we need to put the positions into a common reference frame.

if keyword_set(chip1) then begin
   xycom=transpose([2048d,-18d])
endif else if keyword_set(chip2) then begin
   xycom=transpose([2048d,2065d])
endif

; Correct the common position
xycom_corr=coord_warp(xycom,$
                      chip1=keyword_set(chip1),chip2=keyword_set(chip2),$
                      f775w=keyword_set(f775w),f475w=keyword_set(f475w))

xy_pos_corr=coord_warp(xy_pos,$
                       chip1=keyword_set(chip1),chip2=keyword_set(chip2),$
                       f775w=keyword_set(f775w),f475w=keyword_set(f475w))

xy_twochip=xy_pos_corr
xy_twochip[*,0]-=xycom_corr[0]
xy_twochip[*,1]-=xycom_corr[1]

;print,'COORD_TWOCHIP',size(xy_twochip,/dimensions)

return,xy_twochip

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function acs_geom_distortion,XY_POS,$
                     CHIP1=CHIP1,CHIP2=CHIP2,$
                     F775W=F775W,F475W=F475W

; Check if there is only one position being called
pos_dims=size(xy_pos,/dimensions)
onexy=(n_elements(pos_dims) eq 1) or (pos_dims[0] eq 1)

; Transpose if necessary
dotrans=onexy and (pos_dims[0] ne 1)

if dotrans then xy_pos=transpose(xy_pos)

xy_twochip=coord_twochip(xy_pos,$
                       chip1=keyword_set(chip1),chip2=keyword_set(chip2),$
                       f775w=keyword_set(f775w),f475w=keyword_set(f475w))

xy_twochip_zero=coord_twochip(transpose([0d,0d]),/chip2,$
                              f775w=keyword_set(f775w),$
                              f475w=keyword_set(f475w))

xy_out=xy_pos
xy_out[*,0]=xy_twochip[*,0]-xy_twochip_zero[0,0]
xy_out[*,1]=xy_twochip[*,1]-xy_twochip_zero[0,1]

if dotrans then xy_out=transpose(xy_out)

; De-scale it so we get pixel units
xy_out/=0.05d

return, xy_out
end
