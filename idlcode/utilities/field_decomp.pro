FUNCTION FIELD_DECOMP, DATA_ARRAY, CENTER, $
                       DATAFILE=DATAFILE, XY_COLS=XY_COLS, $
                       FIELD_COLS=FIELD_COLS, SPIN=SPIN, $
                       TO_POLAR=TO_POLAR
  
; This function converts a Cartesian (Fx,Fy) field description into an
; azimuthal (F1,F2) description assuming a given spin symmetry.  The
; output array has the same structure as the input array, except that
; the field array columns are replaced with the projected field
; values. Also, the coordinates are converted into radius,angle format
; if the TO_POLAR keyword is active.

if n_elements(center) ne 2 then center=[0d,0d]

if not keyword_set(spin) then spin=1
spin=long(spin)
if (spin lt 1) then spin=1L
spin=double(spin)

; Input control
if keyword_set(datafile) then data_array=read_data(datafile)
dims=size(data_array,/dimensions)



if not keyword_set(xy_cols) then xy_cols=[0,1]
if not keyword_set(field_cols) then field_cols=[2,3]


npts=dims[0]

dx=data_array[*,xy_cols[0]]-center[0]
dy=data_array[*,xy_cols[1]]-center[1]

fx=data_array[*,field_cols[0]]
fy=data_array[*,field_cols[1]]

if (xy_cols[0] eq field_cols[0]) and (xy_cols[1] eq field_cols[1]) then begin
   fx=dx
   fy=dy
endif

; polar position
radius=sqrt(dx^2+dy^2)
angle=atan(dy,dx)

; polar field
fmag=sqrt(fx^2+fy^2)
fang=atan(fy,fx)

; project onto the spin basis
f1=fmag*cos(fang-spin*angle)
f2=fmag*sin(fang-spin*angle)

out_array=data_array

if keyword_set(to_polar) then begin
   out_array[*,xy_cols[0]]=radius
   out_array[*,xy_cols[1]]=angle
endif

out_array[*,field_cols[0]]=f1
out_array[*,field_cols[1]]=f2

return, out_array

end
