pro compare_pars, p_test, p_ref, tol=tol, ae=ae, pol=pol

; This function compares a test set of parameters to a reference set
; to tell how different they are.  We normalize I0, A/alpha, eps, and
; xi to the reference parameter values.  For each of the fields with
; spin dependence (Xc-Yc, e+-ex, g, G1, G3) we normalize by
; the sum in quadrature of the field.  If the reference value is less
; than the given tolerance, the absolute difference is given.

if not keyword_set(tol) then tol=1d-7

if not (keyword_set(ae) or keyword_set(pol)) then begin
   print,'Choose parameter set...'
   return
endif

diff=p_test-p_ref

ref=p_ref

ref[1:2]=sqrt(total(p_ref[1:2]^2))
ref[6:7]=sqrt(total(p_ref[6:7]^2))
ref[8:9]=sqrt(total(p_ref[8:9]^2))
ref[10:11]=sqrt(total(p_ref[10:11]^2))

if (keyword_set(new) or keyword_set(pol)) then $
   ref[4:5]=sqrt(total(p_ref[4:5]^2))

if keyword_set(ae) then labels=$
   ['I0   ','Xc   ','Yc   ','A    ','eps  ','xi   ',$
    'g1   ','g2   ','G11  ','G12  ','G31  ','G32  ']
if keyword_set(pol) then labels=$
   ['I0   ','Xc   ','Yc   ','A    ','e+   ','ex   ',$
    'g1   ','g2   ','G11','G12','G31','G32']

tolok=where(ref gt tol)

out=strarr(5,14)
out[0,0]='Par  '
out[0,1]='-----'
out[1,0]='  Test     '
out[1,1]='-----------'
out[2,0]='   Ref     '
out[2,1]='-----------'
out[3,0]='  Diff     '
out[3,1]='-----------'
out[4,0]='abs/rel'
out[4,1]='-------'

out[0,2:13]=labels
out[4,2:13]='  abs  '
out[4,tolok+2]='  rel  '

diff[tolok]/=ref[tolok]

out[1,2:13]=string(p_test,format='(e11.3)')
out[2,2:13]=string(p_ref,format='(e11.3)')
out[3,2:13]=string(diff,format='(e11.3)')

print,out

end
