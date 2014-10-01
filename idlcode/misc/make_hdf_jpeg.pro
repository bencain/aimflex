;+
; NAME:
;   make_hdf_jpeg
; BUGS:
;   everything hard wired
;-
pro make_hdf_jpeg

; set parameters
scales= [7000.,4000.,15000.]
nonlinearity= 3.0
resizefactor= 0.5

; read data
im= mrdfits('f814_mosaic_*.fits')
tmp= size(im,/dimensions)
nx= tmp[0]
ny= tmp[1]
RGBim= fltarr(nx,ny,3)
RGBim[*,*,0]= im
RGBim[*,*,1]= mrdfits('f606_mosaic_*.fits')
RGBim[*,*,2]= mrdfits('f450_mosaic_*.fits')

; rebin
RGBim= rebin(RGBim,floor(nx*resizefactor),(ny*resizefactor),3)

; scale and set colors
RGBim = nw_scale_rgb(RGBim,scales=scales)
RGBim = nw_arcsinh_fit(RGBim,nonlinearity=nonlinearity)
RGBim = nw_fit_to_box(RGBim,origin=origin)
RGBim = nw_float_to_byte(RGBim)

; write
WRITE_JPEG,'hdf_wherry.jpg',RGBim,TRUE=3,QUALITY=100
return
end
