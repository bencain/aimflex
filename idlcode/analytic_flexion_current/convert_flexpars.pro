function convert_flexpars, pars, gn_to_psin=gn_to_psin, psin_to_gn=psin_to_gn

; This function converts the flexion fields between the G1/G3
; parametrization and the psi_1/psi_3 parametrization, both as defined
; by Schneider & Er (2008).

outpars=pars
npar=n_elements(pars)

if keyword_set(gn_to_psin) then begin

   g=dcomplex(pars[npar-6],pars[npar-5])
   G1=dcomplex(pars[npar-4],pars[npar-3])
   G3=dcomplex(pars[npar-2],pars[npar-1])

   psi1=(G1-g*conj(G1))/(4d*(1d - conj(g)*g))
   psi3=G3/4d - g*psi1

   outpars[npar-4]=real_part(psi1)
   outpars[npar-3]=imaginary(psi1)
   outpars[npar-2]=real_part(psi3)
   outpars[npar-1]=imaginary(psi3)
endif

if keyword_set(psin_to_gn) then begin

   g=dcomplex(pars[npar-6],pars[npar-5])
   psi1=dcomplex(pars[npar-4],pars[npar-3])
   psi3=dcomplex(pars[npar-2],pars[npar-1])

   G1=4d*(psi1+g*conj(psi1))
   G3=4d*(psi3+g*psi1)

   outpars[npar-4]=real_part(G1)
   outpars[npar-3]=imaginary(G1)
   outpars[npar-2]=real_part(G3)
   outpars[npar-1]=imaginary(G3)

endif

return,outpars

end
