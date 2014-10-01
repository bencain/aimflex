function bnucgs,lambda,temp

@constants.include

; Calculate the Planck function Bnu = 2h nu^3/(exp(h nu/kT)-1)/c^2
; in cgs units
;
; lambda in cm, temp in K
; bnu returned in erg/s/cm/cm/Hz/Str
 

  nu = clight / (lambda ) ; convert to frequency


  c1 = 2*hplanck/clight^2 	; first radiation constant
  c2 = hplanck/kboltz		; second radiation constant

  x = c2 * nu / temp
  bn =  c1*nu^3 / ( exp(x < 87.0) - 1. )

;  print,nu, c1, c2, x

  return, bn


;  if n_elements(temp eq 1) then begin
;
;     bn = lambda* 0.0
;  endif else begin






;  good = where( x lt 88.0, ngood)

;  if( ngood gt 0) then  begin
;
;      bn(good) =  c1*nu(good)^3 / ( exp(x(good)) - 1. )
;
;
;     if n_elements(nu gt 1) then begin 


;     endif else begin
    
;        bn[good] =  c1*nu^3 / ( exp(x[good]) - 1. )
;     endelse


;   endif 


 
 
  end
