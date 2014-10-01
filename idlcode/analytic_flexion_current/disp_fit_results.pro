PRO DISP_FIT_RESULTS, filetag, fitcat, wnum, $
                      disp_obj=disp_obj, runthrough=runthrough, logscale=logscale

; This procedure will display the data image, the fit image, the
; residual image and the desired field over the top of the data image.

; Grab a big window
if ((n_elements(wnum) eq 0) or  $
    (!d.window eq -1) or $
    (set_intersection(!d.window,wnum) eq -1)) then begin
   window,xsize=750,ysize=750, title=title,/free
   wnum=!d.window
endif else begin
   wset,wnum
endelse

ncat=(size(fitcat,/dimensions))[0]
if n_elements(disp_obj) eq 0 then disp_obj=lindgen(ncat)
max_count=floor(alog10(max(fitcat[*,0])))+1

symsize=(0.02d)*375d

; Different spin symbols
;Spin-1: Arrow
s1sym_x=[-1d,1d,(1d)-sqrt(2d)/3d,1d,(1d)-sqrt(2d)/3d]
s1sym_y=[0d,0d,sqrt(2d)/3d,0d,-sqrt(2d)/3d]

;Spin-2: Bar
s2sym_x=[-1d,1d]
s2sym_y=[0d,0d]

;Spin-3: 3-prong
s3sym_x=[1d,0d,-0.5d,0d,-0.5d]
s3sym_y=[0d,0d,sqrt(3d)/2d,0d,-sqrt(3d)/2d]


npar=(size(fitcat,/dimensions))[1]-3

for i=0,n_elements(disp_obj)-1 do begin
   cur_count=floor(alog10(fitcat[disp_obj[i],0]))+1
   gap=''
   for j=cur_count,max_count-1 do gap+='0'
   tail=strcompress(string(long(fitcat[disp_obj[i],0])),/remove_all)+'.fits'
; Read in the images
   dimg=mrdfits(filetag+'_dimg'+gap+tail)
   fimg=mrdfits(filetag+'_fimg'+gap+tail)
   rimg=mrdfits(filetag+'_rimg'+gap+tail)

   goodpix=where(fimg gt 0d,ngoodpix)
   if ngoodpix gt 0 then minval=min(fimg[goodpix]) else minval=0d

   dimg[goodpix]-=minval
   fimg[goodpix]-=minval

   if keyword_set(logscale) then begin
      dimg=alog(dimg+1d)
      fimg=alog(fimg+1d)
   endif

   pars=fitcat[disp_obj[i],3:npar+2]

   imsize=size(dimg,/dimensions)
   scaling=375d/max(imsize)

   tvscl,scale_image(dimg,scaling),0
   tvscl,scale_image(fimg,scaling),1
   tvscl,scale_image(rimg,scaling),2
   tvscl,scale_image(dimg,scaling),3

   !p.multi=[4,2,2]
   ctr_x=fitcat[disp_obj[i],4]
   ctr_y=fitcat[disp_obj[i],5]
   
   gmag=sqrt(total(pars[npar-6:npar-5]^2))
   gphi=atan(pars[npar-5],pars[npar-6])/2d

   psi1mag=sqrt(total(pars[npar-4:npar-3]^2))/scaling
   psi1phi=atan(pars[npar-3],pars[npar-4])/2d

   psi3mag=sqrt(total(pars[npar-2:npar-1]^2))/scaling
   psi3phi=atan(pars[npar-1],pars[npar-2])/2d

   sym1_x=(cos(psi1phi[0])*s1sym_x-sin(psi1phi[0])*s1sym_y)*psi1mag[0]*imsize[0]/0.02d
   sym1_y=(cos(psi1phi[0])*s1sym_y+sin(psi1phi[0])*s1sym_x)*psi1mag[0]*imsize[0]/0.02d

   sym2_x=(cos(gphi[0])*s2sym_x-sin(gphi[0])*s2sym_y)*gmag[0]*imsize[0]/6d
   sym2_y=(cos(gphi[0])*s2sym_y+sin(gphi[0])*s2sym_x)*gmag[0]*imsize[0]/6d

   sym3_x=(cos(psi3phi[0])*s3sym_x-sin(psi3phi[0])*s3sym_y)*psi3mag[0]*imsize[0]/0.02d
   sym3_y=(cos(psi3phi[0])*s3sym_y+sin(psi3phi[0])*s3sym_x)*psi3mag[0]*imsize[0]/0.02d

   plot,sym1_x+ctr_x,sym1_y+ctr_y,$
        xrange=[-1d,1d]*(imsize[0]-1d)/2d,xstyle=5,$
        yrange=[-1d,1d]*(imsize[1]-1d)/2d,ystyle=5,$
        color=fsc_color('blue'), thick=3d
   oplot,sym2_x+ctr_x,sym2_y+ctr_y,$
         color=fsc_color('red'), thick=3d
   oplot,sym3_x+ctr_x,sym3_y+ctr_y,$
         color=fsc_color('orange'), thick=3d
   xyouts,0.05,0.54,$
          'Catalog Number: '+strcompress(string(long(fitcat[disp_obj[i],0])),/remove_all),$
          /normal,charsize=2d,color=fsc_color('green')
   xyouts,0.05,0.51,$
          'Located at (x,y) = ( '+$
          strcompress(string(long(fitcat[disp_obj[i],1])),/remove_all)+', '+$
          strcompress(string(long(fitcat[disp_obj[i],2])),/remove_all)+$
          ' )',/normal,charsize=2d,color=fsc_color('green')
   xyouts,0.55,0.51,'Fit Image',charsize=2d,color=fsc_color('red'),/normal
   xyouts,0.05,0.01,'Residual Image',charsize=2d,color=fsc_color('red'),/normal
   xyouts,0.55,0.01,'Data Image',charsize=2d,color=fsc_color('red'),/normal
   
   if keyword_set(runthrough) then begin 
      wait,1.5
      goon=1
   endif else goon=get_num_resp('Next image? (1=yes, 0=quit): ')
   if goon eq 0 then break
endfor

!p.multi=0
;clear_windows
end
