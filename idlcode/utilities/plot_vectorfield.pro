pro plot_vectorfield, data_array, datafile=datafile, psfile=psfile, $
                      xy_cols=xy_cols, vec_cols=vec_cols,title=title, $
                      select_col=select_col, limit=limit

; This procedure plots a vectorfield from a data array.  If the
; DATAFILE keyword is set, data_array input is ignored and data is
; read in from the file name given.  If the PSFILE keyword is set,
; then the plot is also written to the given postscript file.  XY_COLS
; and VEC_COLS are 2 element arrays indicating where the position and
; vector components are located in data_array,
; e.g. X=data_array[*,xy_cols[0]], etc.  If these are not set, they
; are assumed to be XY_COLS=[0,1] and VEC_COLS=[2,3]

if keyword_set(datafile) then data_array=read_data(datafile)
if not keyword_set(xy_cols) then xy_cols=[0,1]
if not keyword_set(vec_cols) then vec_cols=[2,3]
if not keyword_set(title) then title='Vector Plot'

npts=(size(data_array,/dimensions))[0]

if (keyword_set(select_col) and keyword_set(limit))then begin
   good=where(data_array[*,select_col] lt limit,npts)
endif else good=dindgen(npts)

window,xsize=600,ysize=600, title=title

; Read in the data
x_arr=data_array[good,xy_cols[0]]
y_arr=data_array[good,xy_cols[1]]
vec_x_arr=data_array[good,vec_cols[0]]
vec_y_arr=data_array[good,vec_cols[1]]

xmax=max(x_arr)
ymax=max(y_arr)
xmin=min(x_arr)
ymin=min(y_arr)

xwidth=xmax-xmin
xctr=(xmax+xmin)/2d
ywidth=ymax-ymin
yctr=(ymax+ymin)/2d

width=(1.1d)*max([xwidth,ywidth])
xrange=xctr+width*[-0.5d,0.5d]
yrange=yctr+width*[-0.5d,0.5d]

ignore_val=max(abs([vec_x_arr,vec_y_arr]))+10d

; Scale the peak vector size so that we can see the arrows
vecmags=sqrt(vec_x_arr^2+vec_y_arr^2)
maxlen=(0.01d)*max(vecmags)/median(vecmags)

; Plot the field
partvelvec, vec_x_arr, vec_y_arr, x_arr, y_arr, /isotropic, $
            xrange=xrange, yrange=yrange, length=maxlen

; Save the image if PSFILE is set.

if keyword_set(psfile) then begin
   set_plot,'PS'
   device,filename=psfile
   partvelvec, vec_x_arr, vec_y_arr, x_arr, y_arr, /isotropic,$
               xrange=xrange, yrange=yrange, length=maxlen
   device,/close
   set_plot,'X'
endif


end
