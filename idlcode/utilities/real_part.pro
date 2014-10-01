function real_part, x

; This is a kludge because GDL doesn't carry the REAL_PART
; function for some bizarre reason.  It does however have IMAGINARY,
; so we'll use that.
type=size(x,/type)

if (type eq 5) or (type eq 9) then begin
   out=imaginary(dcomplex(0,1)*x)
endif else begin
   out=imaginary(complex(0,1)*x)
endelse

return,out

end
