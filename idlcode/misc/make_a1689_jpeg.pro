;+
; NAME:
;   make_a1689_jpeg
; BUGS:
;   everything hard wired
;-
pro make_a1689_jpeg

; set parameters
scales= [600.,700.,800.]
nonlinearity= 1.0
resizefactor= 0.5

; read data
dir='data/a1689_hla/'
im= mrdfits(dir+'F775W.fits',1)
tmp= size(im,/dimensions)
nx= tmp[0]
ny= tmp[1]
RGBim= fltarr(nx,ny,3)
RGBim[*,*,0]= im
RGBim[*,*,1]= mrdfits(dir+'F625W.fits',1)
RGBim[*,*,2]= mrdfits(dir+'F475W.fits',1)

RGBim=RGBim[1450:4349,1450:4349,*]
tmp=size(RGBim,/dimensions)
nx=tmp[0]
ny=tmp[1]

; rebin
RGBim= rebin(RGBim,floor(nx*resizefactor),(ny*resizefactor),3)

; scale and set colors
RGBim = nw_scale_rgb(RGBim,scales=scales)
RGBim = nw_arcsinh_fit(RGBim,nonlinearity=nonlinearity)
RGBim = nw_fit_to_box(RGBim,origin=origin)
RGBim = nw_float_to_byte(RGBim)

; write
WRITE_JPEG,'a1689_1.jpg',RGBim,TRUE=3,QUALITY=100
return
end
