PRO RADIAL_FLEXVAR, DATA_ARR,CENTER, NBIN=NBIN, NWEDGE=NWEDGE, $
                    XY_COLS=XY_COLS, FLEX_COLS=FLEX_COLS, MAXFLEX=MAXFLEX, $
                    SELECT_COL=SELECT_COL, LIMIT=LIMIT, OK=OK, NOK=NOK,$
                    TO_EPS=TO_EPS, PSTAG=PSTAG, $
                    RAVES=RAVES,TAVES=TAVES, RVARS=RVARS, TVARS=TVARS,$ 
                    NOBJ=NOBJ,WEDGE_PHASE=WEDGE_PHASE, NDENSITY=NDENSITY,$
                    FLEXIONSCALE=FLEXIONSCALE


datasize=size(data_arr,/dimensions)

nok=datasize[0]
ok=lindgen(nok)

if not keyword_set(xy_cols) then xy_cols=[1,2]
if not keyword_set(flex_cols) then flex_cols=[11,12]
if not keyword_set(flexionscale) then flexionscale=1d3
if not keyword_set(wedge_phase) then wedge_phase=0d 
wedge_phase+=!dpi

; Select with a flexion maximum
fmag=sqrt(total(data_arr[*,flex_cols[0]:flex_cols[1]]^2,2))

if not keyword_set(maxflex) then maxflex=1d*flexionscale
ok=set_intersection(ok,where(fmag lt maxflex),count=nok)

; Select with a FOM maximum
if (keyword_set(select_col) and keyword_set(limit)) then $
   ok=set_intersection(ok,where(data_arr[ok,select_col] lt limit),count=nok)

; Make sure there is a center
if n_elements(center) ne 2 then center=[median(data_arr[*,xy_cols[0]]),$
                                        median(data_arr[*,xy_cols[1]])]

xmax=ceil(max(data_arr[*,xy_cols[0]]))
ymax=ceil(max(data_arr[*,xy_cols[1]]))

; Set the number of bins
if not keyword_set(nbin) then nbin=10
if not keyword_set(nwedge) then nwedge=1

; Only go if there are good objects
if nok gt 0 then begin

; Put the flexion data into a radial format.
   proj_data=field_decomp(data_arr[ok,*],center,xy_cols=xy_cols,$
                          field_cols=flex_cols,spin=1, /to_polar)

   radii=proj_data[*,xy_cols[0]]
   angles=proj_data[*,xy_cols[1]]

   rmin=min(radii)
   rmin*=0.95d

;   rmax=max(radii)
   absdx=abs(data_arr[*,xy_cols[0]]-center[0])
   absdy=abs(data_arr[*,xy_cols[1]]-center[1])
   rmax=max([absdx,absdy])
   rmax*=1.05d

   dr=(rmax-rmin)/double(nbin)
   dt=2d*!dpi/nwedge

   rbins=rmin+dr*(dindgen(nbin)+0.5d)

   nobj=dblarr(nbin,nwedge+1)
   ndensity=dblarr(nbin,nwedge+1)
   raves=dblarr(nbin,nwedge+1)
   taves=dblarr(nbin,nwedge+1)
   rvars=dblarr(nbin,nwedge+1)
   tvars=dblarr(nbin,nwedge+1)

   nobj[*,0]=rbins
   ndensity[*,0]=rbins
   raves[*,0]=rbins
   taves[*,0]=rbins
   rvars[*,0]=rbins
   tvars[*,0]=rbins

   nx=100
   ny=100
   img_rad=dblarr(nx,ny)-1d
   img_tan=dblarr(nx,ny)-1d

   ximg=xmax*(dindgen(nx,ny) mod nx)/(nx-1)-center[0]
   yimg=ymax*double(lindgen(nx,ny)/long(ny))/(ny-1)-center[1]
   rimg=sqrt(ximg^2+yimg^2)
   timg=atan(yimg,ximg)

   angles=atan(sin(angles-wedge_phase),cos(angles-wedge_phase))
   timg=atan(sin(timg-wedge_phase),cos(timg-wedge_phase))

   for i=0,nbin-1 do begin
      rlower=rbins[i]-dr/2d
      rupper=rbins[i]+dr/2d
   
      for j=0,nwedge-1 do begin
         tlower=dt*double(j)-!dpi
         tupper=dt*double(j+1)-!dpi
         ruse=where((radii ge rlower) and $
                     (radii lt rupper),nr)
         tuse=where((angles ge tlower) and $
                    (angles le tupper),nt)
         
         use=set_intersection(ruse,tuse,count=nuse)
;         print,nr
;         print,nt
;         print,nuse

         imuse=where((rimg ge rlower) and $
                     (rimg le rupper) and $
                     (timg ge tlower) and $
                     (timg le tupper))

         if nuse gt 3 then begin
            raves[i,j+1]=mean(proj_data[use,flex_cols[0]])
            taves[i,j+1]=mean(proj_data[use,flex_cols[1]])
            rvars[i,j+1]=variance(proj_data[use,flex_cols[0]])
            tvars[i,j+1]=variance(proj_data[use,flex_cols[1]])
         endif
         nobj[i,j+1]=nuse
         ndensity[i,j+1]=double(nuse)/(0.5d*dt*(rupper^2-rlower^2))
         img_rad[imuse]=rvars[i,j+1]
         img_tan[imuse]=tvars[i,j+1]
      endfor
   endfor


   if nwedge gt 1 then begin
      rvmax=max(rvars[*,1:nwedge])
      rvmin=min(rvars[*,1:nwedge]) 
      tvmax=max(tvars[*,1:nwedge])
      tvmin=min(tvars[*,1:nwedge]) 
   endif else begin
      rvmax=max(rvars[*,1])
      rvmin=min(rvars[*,1])
      tvmax=max(tvars[*,1])
      tvmin=min(tvars[*,1])
   endelse

   rvc=(rvmax+rvmin)/2d
   drv=rvmax-rvmin
   tvc=(tvmax+tvmin)/2d
   dtv=tvmax-tvmin
   if drv eq 0d then drv=0.1d
   if dtv eq 0d then dtv=0.1d

   rc=(rmax+rmin)/2d
   
   syms=-[1,2,4,5,6,7]
   lines=[0,2,4,5]
   nsym=n_elements(syms)
   nline=n_elements(lines)

   if keyword_set(pstag) then begin

; Plot the radial (real) component to ps file
      setup_ps,pstag+'_r_flexvar.ps'
      t=5.0d
      plot,rbins,rvars[*,1],xrange=[rc-0.55d*nbin*dr,rc+0.55d*nbin*dr],xstyle=1,$
           yrange=[rvc-0.55d*drv,rvc+0.55d*drv],ystyle=1, linestyle=lines[0],$
           psym=syms[0],xmargin=[17,10], ymargin=[8,8],$
           charsize=2.5,charthick=t,thick=t,xthick=t,ythick=t,$
           xtitle='Radius',ytitle='Radial 1-Flexion Variance (1/pixel)^2'

      for i=1,nwedge-1 do begin
         oplot,rbins,rvars[*,i+1],linestyle=lines[i mod nline],$
               psym=syms[i mod nsym],thick=t
      endfor

      setup_x

; Plot the tangential (imaginary) component to ps file
      setup_ps,pstag+'_t_flexvar.ps'
      t=5.0d
      plot,rbins,tvars[*,1],xrange=[rc-0.55d*nbin*dr,rc+0.55d*nbin*dr],xstyle=1,$
           yrange=[tvc-0.55d*dtv,tvc+0.55d*dtv],ystyle=1, linestyle=lines[0],$
           psym=syms[0],xmargin=[17,10], ymargin=[8,8],$
           charsize=2.5,charthick=t,thick=t,xthick=t,ythick=t,$
           xtitle='Radius',ytitle='Tangential 1-Flexion Variance (1/pixel)^2'

      for i=1,nwedge-1 do begin
         oplot,rbins,tvars[*,i+1],linestyle=lines[i mod nline],$
               psym=syms[i mod nsym],thick=t
      endfor

      setup_x


      if keyword_set(to_eps) then begin
         ps2eps,pstag+'_r_flexvar.ps',/delete
         ps2eps,pstag+'_t_flexvar.ps',/delete
      endif

      mwrfits,img_rad,pstag+'_r_flexvar_img.fits'
      mwrfits,img_tan,pstag+'_t_flexvar_img.fits'

   endif else begin
      
;Plot the real/radial component to window
      window,/free,xsize=500,ysize=800
      !p.multi=[0,1,2]
      plot,rbins,rvars[*,1],xrange=[rc-0.55d*nbin*dr,rc+0.55d*nbin*dr],$
           xstyle=1,yrange=[rvc-0.55d*drv,rvc+0.55d*drv],ystyle=1, $
           linestyle=lines[0],psym=syms[0],xmargin=[17,10], ymargin=[8,8],$
           charsize=2.0,charthick=2.0,thick=2.0,xthick=2.0,ythick=2.0,$
           xtitle='Radius',ytitle='Radial 1-Flexion Variance (1/pixel)^2'

      for i=1,nwedge-1 do begin
         oplot,rbins,rvars[*,i+1],linestyle=lines[i mod nline],$
               psym=syms[i mod nsym],thick=2.0
      endfor

;Plot the imaginary/tangential component to window
      plot,rbins,tvars[*,1],xrange=[rc-0.55d*nbin*dr,rc+0.55d*nbin*dr],$
           xstyle=1,yrange=[tvc-0.55d*dtv,tvc+0.55d*dtv],ystyle=1, $
           linestyle=lines[0],psym=syms[0],xmargin=[17,10], ymargin=[8,8],$
           charsize=2.0,charthick=2.0,thick=2.0,xthick=2.0,ythick=2.0,$
           xtitle='Radius',ytitle='Tangential 1-Flexion Variance (1/pixel)^2'

      for i=1,nwedge-1 do begin
         oplot,rbins,tvars[*,i+1],linestyle=lines[i mod nline],$
               psym=syms[i mod nsym],thick=2.0
      endfor

      !p.multi=0

; Plot the number/number density of objects
      window,/free,xsize=500,ysize=800
      !p.multi=[0,1,2]
      
      if nwedge le 1 then nmax=max(nobj[*,1]) else $
         nmax=max(nobj[*,1:nwedge-1])
      plot,nobj[*,0],nobj[*,1],xrange=[rc-0.55d*nbin*dr,rc+0.55d*nbin*dr],$
           xstyle=1, yrange=[0d,1.1d*nmax],ystyle=1,$
           linestyle=lines[0],psym=syms[0],xmargin=[17,10], ymargin=[8,8],$
           charsize=2.0,charthick=2.0,thick=2.0,xthick=2.0,ythick=2.0,$
           xtitle='Radius',ytitle='Object Number'

      for i=1,nwedge-1 do begin
         oplot,rbins,nobj[*,i+1],linestyle=lines[i mod nline],$
               psym=syms[i mod nsym],thick=2.0
      endfor

      if nwedge le 1 then ndmax=max(ndensity[*,1]) else $
         ndmax=max(ndensity[*,1:nwedge-1])
      plot,ndensity[*,0],ndensity[*,1],$
           xrange=[rc-0.55d*nbin*dr,rc+0.55d*nbin*dr],$
           xstyle=1, yrange=[0d,1.1d*ndmax],ystyle=1,$
           linestyle=lines[0],psym=syms[0],xmargin=[17,10], ymargin=[8,8],$
           charsize=2.0,charthick=2.0,thick=2.0,xthick=2.0,ythick=2.0,$
           xtitle='Radius',ytitle='Object Number Density'

      for i=1,nwedge-1 do begin
         oplot,rbins,ndensity[*,i+1],linestyle=lines[i mod nline],$
               psym=syms[i mod nsym],thick=2.0
      endfor

      !p.multi=0

      disp_scaled_image,img_rad,wnum_rad,title='Radial 1-Flexion Variance'
      disp_scaled_image,img_tan,wnum_tan,title='Tangential 1-Flexion Variance'
   endelse


endif else print,'No usable objects: all above MAXFLEX...'

end
