FUNCTION POWERLAW_LENSMODEL, CMPLX_POS, A, K, INCLUDE_KAPPA=INCLUDE_KAPPA,$
                             FLEXIONSCALE=FLEXIONSCALE

; This takes in a complex position (relative to the lens center), the
; amplitude of the lens potential and the power of the potental.  We
; assume an azimuthally symmetric lens potential of the form.
;
;    psi= A * (cmplx_pos * conj(cmplx_pos))^K
;

if not keyword_set(flexionscale) then flexionscale=1d3

kappa=2d*A*K^2*(cmplx_pos*conj(cmplx_pos))^(K-1d)
gamma=2d*A*K*(K-1d)*cmplx_pos^2*(cmplx_pos*conj(cmplx_pos))^(K-2d)

fflex=4d*A*K^2*(K-1d)*cmplx_pos*(cmplx_pos*conj(cmplx_pos))^(K-2d)
gflex=4d*A*K*(K-1d)*(K-2d)*cmplx_pos^3*(cmplx_pos*conj(cmplx_pos))^(K-3d)

n=min([n_elements(kappa),$
       n_elements(gamma),$
       n_elements(fflex),$
       n_elements(gflex)])

shear=gamma/(1d -kappa)
G1=(fflex+shear*conj(fflex))/(1d -kappa)
G3=(gflex+shear*fflex)/(1d -kappa)

lensmodel=dblarr(n,8)

lensmodel[*,0]=real_part(cmplx_pos)
lensmodel[*,1]=imaginary(cmplx_pos)

lensmodel[*,2]=real_part(shear)
lensmodel[*,3]=imaginary(shear)

lensmodel[*,4]=real_part(G1)*flexionscale
lensmodel[*,5]=imaginary(G1)*flexionscale

lensmodel[*,6]=real_part(G3)*flexionscale
lensmodel[*,7]=imaginary(G3)*flexionscale

if keyword_set(include_kappa) then begin
   out=dblarr(n,9)
   out[*,0:7]=lensmodel
   out[*,8]=kappa
endif else out=lensmodel

if n lt 2 then return,transpose(out) else return,out

END
