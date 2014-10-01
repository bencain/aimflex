pro plot_field, data_array, wnum, $
                datafile=datafile, xy_cols=xy_cols, spin=spin, $
                field_cols=field_cols, select_col=select_col, limit=limit,$
                symscale=symscale, log=log, psfile=psfile, title=title,$
                autoscale=autoscale, overplot=overplot,$
                standard=standard,std_x=std_x,std_y=std_y,std_label=std_label

; This is an update to PLOT_VECTORFIELD which plots a 2-D complex
; field.  The SPIN keyword determines whether the field is plotted as
; a vector, a polarization or as a spin-3 field.  If PSFILE is set,
; then the plot is also written to that postscript file instead of to
; screen. 

if not keyword_set(xy_cols) then xy_cols=[0,1]
if not keyword_set(field_cols) then field_cols=[2,3]
if not keyword_set(title) then title='Field Plot'

if not keyword_set(psfile) then begin
   if ((n_elements(wnum) eq 0) or  (!d.window eq -1)) then begin
      window,xsize=750,ysize=750, title=title,/free
      wnum=!d.window
   endif else begin
      wset,wnum
   endelse
endif

if not keyword_set(spin) then spin=1
spin=long(spin)
if (spin lt 1) or (spin gt 3) then spin=1L

; Input control
if keyword_set(datafile) then data_array=read_data(datafile)
dims=size(data_array,/dimensions)
if dims[1] lt 4 then begin
   fillarr=dblarr(dims[0],4)
   fillarr[*,0:dims[1]-1]=data_array
   data_array=fillarr
endif

npts=dims[0]

if (keyword_set(select_col) and keyword_set(limit))then begin
   good=where(data_array[*,select_col] lt limit,npts)
endif else good=dindgen(npts)

if not keyword_set(symscale) then symscale=1d

; Read in the data
x_arr=data_array[good,xy_cols[0]]
y_arr=data_array[good,xy_cols[1]]
field_x_arr=data_array[good,field_cols[0]]
field_y_arr=data_array[good,field_cols[1]]

xmax=max(x_arr)
ymax=max(y_arr)
xmin=min(x_arr)
ymin=min(y_arr)

; Find the XY scale of the plot
xwidth=xmax-xmin
xctr=(xmax+xmin)/2d
ywidth=ymax-ymin
yctr=(ymax+ymin)/2d

;xwidth=max([xwidth,30d])
;ywidth=max([ywidth,30d])

width=max([xwidth,ywidth]);*1.1d
xrange=xctr+width*[-0.5d,0.5d]
yrange=yctr+width*[-0.5d,0.5d]

; Set the basic scale of the symbols with respect to the plot.
symsize=(0.02d)*symscale*width

; Create the basic symbols for each of the spin cases

;Spin-1: Arrow
s1sym_x=symsize*[-1d,1d,(1d)-sqrt(2d)/3d,1d,(1d)-sqrt(2d)/3d]
s1sym_y=symsize*[0d,0d,sqrt(2d)/3d,0d,-sqrt(2d)/3d]

;Spin-2: Bar
s2sym_x=symsize*[-1d,1d]
s2sym_y=symsize*[0d,0d]

;Spin-3: 3-prong
s3sym_x=symsize*[1d,0d,-0.5d,0d,-0.5d]
s3sym_y=symsize*[0d,0d,sqrt(3d)/2d,0d,-sqrt(3d)/2d]


;Now plot the field

fieldmag=sqrt(field_x_arr^2+field_y_arr^2)
minfield=min(fieldmag)
if minfield eq 0d then minfield=1d-8
if keyword_set(log) then fieldmag=alog(fieldmag/minfield)

if keyword_set(autoscale) then fieldmag/=mean(fieldmag)

fieldphi=atan(field_y_arr,field_x_arr)/double(spin)

case spin of
   1: begin
      sym_x=(cos(fieldphi[0])*s1sym_x-sin(fieldphi[0])*s1sym_y)*fieldmag[0]
      sym_y=(cos(fieldphi[0])*s1sym_y+sin(fieldphi[0])*s1sym_x)*fieldmag[0]
   end
   2: begin
      sym_x=(cos(fieldphi[0])*s2sym_x-sin(fieldphi[0])*s2sym_y)*fieldmag[0]
      sym_y=(cos(fieldphi[0])*s2sym_y+sin(fieldphi[0])*s2sym_x)*fieldmag[0]
   end
   3: begin
      sym_x=(cos(fieldphi[0])*s3sym_x-sin(fieldphi[0])*s3sym_y)*fieldmag[0]
      sym_y=(cos(fieldphi[0])*s3sym_y+sin(fieldphi[0])*s3sym_x)*fieldmag[0]
   end
endcase

if keyword_set(psfile) then begin
   setup_ps,psfile
   t=5
   ct=5
endif else begin
   t=2
   ct=1
endelse

if keyword_set(overplot) then begin
   oplot,x_arr[0]+sym_x,y_arr[0]+sym_y,thick=t
endif else begin
   plot,x_arr[0]+sym_x,y_arr[0]+sym_y,$
     xrange=xrange,yrange=yrange,/isotropic,thick=t,charsize=2,charthick=ct,xstyle=1,ystyle=1
endelse

for i=1,npts-1 do begin

   case spin of
      1: begin
         sym_x=$
            (cos(fieldphi[i])*s1sym_x-sin(fieldphi[i])*s1sym_y)*fieldmag[i]
         sym_y=$
            (cos(fieldphi[i])*s1sym_y+sin(fieldphi[i])*s1sym_x)*fieldmag[i]
      end
      2: begin
         sym_x=$
            (cos(fieldphi[i])*s2sym_x-sin(fieldphi[i])*s2sym_y)*fieldmag[i]
         sym_y=$
            (cos(fieldphi[i])*s2sym_y+sin(fieldphi[i])*s2sym_x)*fieldmag[i]
      end
      3: begin
         sym_x=$
            (cos(fieldphi[i])*s3sym_x-sin(fieldphi[i])*s3sym_y)*fieldmag[i]
         sym_y=$
            (cos(fieldphi[i])*s3sym_y+sin(fieldphi[i])*s3sym_x)*fieldmag[i]
      end
   endcase

   oplot,x_arr[i]+sym_x,y_arr[i]+sym_y,thick=t

endfor

if keyword_set(standard) then begin
; Put in a standard for reference
   if not keyword_set(std_x) then std_x=xrange[0]+0.75d*xwidth
   if not keyword_set(std_y) then std_y=yrange[0]+0.10d*ywidth
   case spin of
      1: begin
         sym_x=s1sym_x
         sym_y=s1sym_y
      end
      2: begin
         sym_x=s2sym_x
         sym_y=s2sym_y
      end
      3: begin
         sym_x=s3sym_x
         sym_y=s3sym_y
      end
   endcase

   loadct,13
   oplot,std_x+sym_x*10d,std_y+sym_y*10d
   if keyword_set(std_label) then begin
      labelx=xrange[0]+0.1d*xwidth
      xyouts,labelx,std_y,'Field = 10'
   endif
   loadct,0
endif

if keyword_set(psfile) then setup_x

end



         
