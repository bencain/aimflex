function rp_parameter_draw,pmin,pmax,seed,$
                           ctr_lim=ctr_lim,$
                           e_lim=e_lim, g_lim=g_lim,$
                           psi1_lim=psi1_lim, psi3_lim=psi3_lim,$
                           sersic=sersic

if keyword_set(sersic) then npar=13 else npar=12
if not keyword_set(ctr_lim) then ctr_lim=5d
if not keyword_set(e_lim) then e_lim=1d
if not keyword_set(g_lim) then g_lim=1d
if not keyword_set(psi1_lim) then psi1_lim=0.02d
if not keyword_set(psi3_lim) then psi3_lim=0.02d

; This function will take in an array of maxima and minima for a
; rectangular parameter space, plus an optional set of limits for the
; paired parameters and will draw until a set of parameters meets all
; the criteria.  There are default limits if none are input.

newpar=dblarr(npar)

ctrok=(0 eq 1)
eok=(0 eq 1)
gok=(0 eq 1)
psi1ok=(0 eq 1)
psi3ok=(0 eq 1)

allok=(0 eq 1)

; Draw the logN0 and alpha numbers - these have no other limits
newpar[0]=(pmax-pmin)[0]*randomu(seed)+pmin[0]
newpar[3]=(pmax-pmin)[3]*randomu(seed)+pmin[3]
if npar eq 13 then newpar[6]=(pmax-pmin)[6]*randomu(seed)+pmin[6]

; Draw the rest and make sure they're ok.  Redraw if the pair is 
; outside the mag limits.
while not allok do begin
   
   if not ctrok then begin
      ctr=(pmax-pmin)[1:2]*randomu(seed,2)+pmin[1:2]
      newpar[1:2]=ctr
      ctrok=sqrt(total(ctr^2)) lt ctr_lim
   endif

   if not eok then begin
      e=(pmax-pmin)[4:5]*randomu(seed,2)+pmin[4:5]
      newpar[4:5]=e
      eok=sqrt(total(e^2)) lt e_lim
   endif

   if not gok then begin
      g=(pmax-pmin)[npar-6:npar-5]*randomu(seed,2)+pmin[npar-6:npar-5]
      newpar[npar-6:npar-5]=g
      gok=sqrt(total(g^2)) lt g_lim
   endif

   if not psi1ok then begin
      psi1=(pmax-pmin)[npar-4:npar-3]*randomu(seed,2)+pmin[npar-4:npar-3]
      newpar[npar-4:npar-3]=psi1
      psi1ok=sqrt(total(psi1^2)) lt psi1_lim
   endif

   if not psi3ok then begin
      psi3=(pmax-pmin)[npar-2:npar-1]*randomu(seed,2)+pmin[npar-2:npar-1]
      newpar[npar-2:npar-1]=psi3
      psi3ok=sqrt(total(psi3^2)) lt psi3_lim
   endif
   
   allok=ctrok and eok and gok and psi1ok and psi3ok

endwhile

return, newpar

end
