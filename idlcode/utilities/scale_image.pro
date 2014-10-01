FUNCTION SCALE_IMAGE,IMG,SCALE

in_size=size(img)
out_dims=ceil(in_size[1:2]*scale)
in_x=dindgen(in_size[1:2]) mod in_size[1]
in_y=double(lindgen(in_size[1:2])/long(in_size[1]))

out=make_array(out_dims,type=in_size[3])

if scale lt 1d then begin
;   for i=0,out_dims[0]-1 do begin
;      for j=0,out_dims[1]-1 do begin
;         
;         out[i,j]=mean(img[where((i/scale lt in_x+1) and $
;                                 ((i+1)/scale gt in_x) and $
;                                 (j/scale lt in_y+1) and $
;                                 ((j+1)/scale gt in_y))])
;      endfor
;   endfor
   out=img
endif else begin

;;;;;;;;;;;;;;;;;;
; Old Version;;;;;
;;;;;;;;;;;;;;;;;;
;in_size=size(img)
;
;out=make_array(ceil(in_size[1:2]*scale),type=in_size[3])
;
   for i=0,in_size[1]-1 do begin
      for j=0,in_size[2]-1 do begin
         
         out[floor(i*scale):ceil((i+1)*scale)-1,$
             floor(j*scale):ceil((j+1)*scale)-1]=img[i,j]
      endfor
   endfor
endelse

return,out

end
