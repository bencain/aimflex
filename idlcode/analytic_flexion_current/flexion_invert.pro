function flexion_invert, X, PSI1

; This function does a Fourier inversion of the reduced 1-flexion and
; to yeild a mass map. X and PSI1 are complex and assumed to be NxN
; arrays. 

j=dcomplex(0d,1d)

n=size(x,/dimensions)

; Need to find the Fourier ranges.  Dk=2pi/dx and dk=2pi/DX
bigDx1=max(real_part(x))-min(real_part(x))
bigDx2=max(imaginary(x))-min(imaginary(x))

dk1=2d*!dpi/bigDx1
dk2=2d*!dpi/bigDx2

; Make the wavenumbers
k1=shift(dk1*((dindgen(n) mod double(n[0])) - (n[0]-1)/2),0.5d*(n[0]-(n[0] mod 2)) + 1,0)
k2=shift(dk2*(double(lindgen(n)/long(n[0])) - (n[1]-1)/2),0,0.5d*(n[0]-(n[0] mod 2)) + 1)
k=dcomplex(k1,k2)

; Invert!
psi1bar=fft(psi1)

Kbar=(4d*j/k)*psi1bar
kbar[where(abs(k) eq 0)]=0d ; Kill the constant term since we're integrating

bigK=fft(kbar,-1)

kappa=1d - exp(bigK)

; Use the mass-sheet degeneracy to put the min of kappa at zero
lambda=1d/(1d - min(real_part(kappa))) 
kappa_p=lambda*kappa - (1d - lambda)

if min(real_part(kappa_p) ne 0d) then begin
   lambda=1d/(1d - max(real_part(kappa))) 
   kappa_p=lambda*kappa - (1d - lambda)
endif

return,kappa_p

end


