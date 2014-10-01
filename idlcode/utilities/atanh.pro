function atanh,x

; Inverse hyperbolic tangent
;
;
; tanh(y) = sinh(y)/cosh(y) = x
;         = (e^y - e^-y)/(e^y + e^-y)
;         = (e^2y - 1)/(e^2y + 1)
; y = 0.5 alog((1 + x)/(1 - x))
;
; This fails on complex input.

on_error,2
is1=where(abs(abs(x)-1d) lt 1d-10,n1,complement=not1,ncomplement=nnot1)

y=x*0

if n1 gt 0 then begin
   y[is1]=atanh(1d - 1d-9)
endif
if nnot1 gt 0 then begin
   y[not1]=alog((1d + abs(x[not1]))/(1d - abs(x[not1])))/2d
endif

neg=where(x lt 0d, ct)

if ct gt 0 then y[neg]=-y[neg]

return,y

end
