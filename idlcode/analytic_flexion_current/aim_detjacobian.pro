function aim_detjacobian, pars, scale=scale
  
; This function returns the determinant of the lensing Jacobian in the
; quadratic lensing approximation (from Schneider and Er 2008, with
; some further tweaks to put everything into psi_n notation).

npars=n_elements(pars)

if not keyword_set(scale) then scale=pars[3]
g=dcomplex(pars[npars-6],pars[npars-5])
psi1=dcomplex(pars[npars-4],pars[npars-3])
psi3=dcomplex(pars[npars-2],pars[npars-1])

eta=4d*psi1 + 2d*(g*conj(psi1) + conj(g)*psi3)

detA=real_part(1d - g*conj(g)-2d*eta*scale)

return,detA

end
