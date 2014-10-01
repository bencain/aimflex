pro aim_disp_results, dataimage, pars, wnum, $
                     center=center, ellip=ellip, $
                     shear=shear, psi1=psi1, psi3=psi3, $
                     datafile=datafile, psfile=psfile,$
                     fieldscale=fieldscale

; Display the data image and overlay the desired field

if keyword_set(datafile) then dataimage=mrdfits(datafile)

if keyword_set(center) then title='Center Position'
if keyword_set(ellip) then title='Intrinsic Ellipticity'
if keyword_set(shear) then title='Shear'
if keyword_set(psi1) then title='Psi_1'
if keyword_set(psi3) then title='Psi_3'

if keyword_set(psfile) then setup_ps,psfile

disp_scaled_image,dataimage,wnum,title=title
plot,[0,100,100,0],[0,0,100,100],$
     /psym,xstyle=4,ystyle=4,$
     xmargin=[0,0],ymargin=[0,0]
disp_scaled_image,dataimage,wnum

; Symbols
;Spin-1: Arrow
s1sym_x=[-1d,1d,(1d)-sqrt(2d)/3d,1d,(1d)-sqrt(2d)/3d]/2d
s1sym_y=[0d,0d,sqrt(2d)/3d,0d,-sqrt(2d)/3d]/2d

;Spin-2: Bar
s2sym_x=[-1d,1d]/2d
s2sym_y=[0d,0d]

;Spin-3: 3-prong
s3sym_x=[1d,0d,-0.5d,0d,-0.5d]/2d
s3sym_y=[0d,0d,sqrt(3d)/2d,0d,-sqrt(3d)/2d]/2d

npar=n_elements(pars)
ctr=[50d,50d]

if keyword_set(center) then begin

   if not keyword_set(fieldscale) then fieldscale=1d1
   field=pars[1:2]
   phi=atan(field[1],field[0])
   mag=sqrt(total(field^2))
   
   sym_x=(cos(phi)*s1sym_x-sin(phi)*s1sym_y)*mag
   sym_y=(cos(phi)*s1sym_y+sin(phi)*s1sym_x)*mag

   oplot,ctr[0]+sym_x*fieldscale,ctr[1]+sym_y*fieldscale,thick=5d,color=200

endif


if keyword_set(ellip) then begin

   if not keyword_set(fieldscale) then fieldscale=1d2
   field=pars[4:5]
   phi=atan(field[1],field[0])
   mag=sqrt(total(field^2))
   
   sym_x=(cos(phi)*s2sym_x-sin(phi)*s2sym_y)*mag
   sym_y=(cos(phi)*s2sym_y+sin(phi)*s2sym_x)*mag

   oplot,ctr[0]+sym_x*fieldscale,ctr[1]+sym_y*fieldscale,thick=5d,color=200

endif

if keyword_set(shear) then begin

   if not keyword_set(fieldscale) then fieldscale=1d2
   field=pars[npar-6:npar-5]
   phi=atan(field[1],field[0])
   mag=sqrt(total(field^2))
   
   sym_x=(cos(phi)*s2sym_x-sin(phi)*s2sym_y)*mag
   sym_y=(cos(phi)*s2sym_y+sin(phi)*s2sym_x)*mag

   oplot,ctr[0]+sym_x*fieldscale,ctr[1]+sym_y*fieldscale,thick=5d,color=200

endif

if keyword_set(psi1) then begin

   if not keyword_set(fieldscale) then fieldscale=1d4
   field=pars[npar-4:npar-3]
   phi=atan(field[1],field[0])
   mag=sqrt(total(field^2))
   
   sym_x=(cos(phi)*s1sym_x-sin(phi)*s1sym_y)*mag
   sym_y=(cos(phi)*s1sym_y+sin(phi)*s1sym_x)*mag

   oplot,ctr[0]+sym_x*fieldscale,ctr[1]+sym_y*fieldscale,thick=5d,color=200

endif

if keyword_set(psi3) then begin

   if not keyword_set(fieldscale) then fieldscale=1d4
   field=pars[npar-4:npar-3]
   phi=atan(field[1],field[0])
   mag=sqrt(total(field^2))
   
   sym_x=(cos(phi)*s3sym_x-sin(phi)*s3sym_y)*mag
   sym_y=(cos(phi)*s3sym_y+sin(phi)*s3sym_x)*mag

   oplot,ctr[0]+sym_x*fieldscale,ctr[1]+sym_y*fieldscale,thick=5d,color=200

endif

if keyword_set(psfile) then setup_x

end

