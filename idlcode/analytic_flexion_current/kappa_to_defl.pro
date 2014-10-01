FUNCTION KAPPA_TO_DEFL, X, KAPPA

; This helper function converts a covergence map to a deflection
; map. KAPPA can be complex, and DELTA, the deflexion map, is
; complex. X is the array of positions for KAPPA

dx=abs(x[1,0]-x[0,0])
dy=abs(x[0,1]-x[0,0])

s=size(kappa,/dimensions)

; Use the mass sheet degeneracy to make kappa zero-mean
lambda=1d/(1d - mean(kappa))
kappa_p=lambda*kappa + (1d - lambda)

; FFT
kappabar=fft(kappa_p,/double)

kstar=dcomplex((2d*!dpi/(s[0]*dx))*shift((dindgen(s) mod s[0]) - $
                                         0.5d*(s[0]+(s[0] mod 2)) + 1,$
                                         0.5d*(s[0]-(s[0] mod 2)) + 1,0),$
               -(2d*!dpi/(s[1]*dx))*shift(double(lindgen(s)/long(s[0])) - $
                                          0.5d*(s[1]+(s[1] mod 2)) + 1,$
                                          0,0.5d*(s[1]-(s[1] mod 2)) + 1))

deltabar=dcomplex(0d,-2d)*kappabar/kstar
deltabar[0]=0d

delta_p=fft(deltabar,/double,/inverse)

; Now undo the mass sheet transformation
delta=(delta_p + (lambda - 1d)*x)/lambda

return,delta

END

