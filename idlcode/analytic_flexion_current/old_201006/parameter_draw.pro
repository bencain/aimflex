function parameter_draw,pmin,pmax,seed,$
                        ctr_lim=ctr_lim,$
                        e_lim=e_lim, g_lim=g_lim,$
                        G1_lim=G1_lim, G3_lim=G3_lim,$
                        flexionscale=flexionscale

if not keyword_set(flexionscale) then flexionscale=1d3
if not keyword_set(ctr_lim) then ctr_lim=5d
if not keyword_set(e_lim) then e_lim=1d
if not keyword_set(g_lim) then g_lim=1d
if not keyword_set(G1_lim) then G1_lim=0.02d*flexionscale
if not keyword_set(G3_lim) then G3_lim=0.02d*flexionscale

; This function will take in an array of maxima and minima for a
; rectangular parameter space, plus an optional set of limits for the
; paired parameters and will draw until a set of parameters meets all
; the criteria.  There are default limits if none are input.

newpar=dblarr(12)

ctrok=(0 eq 1)
eok=(0 eq 1)
gok=(0 eq 1)
G1ok=(0 eq 1)
G3ok=(0 eq 1)

allok=(0 eq 1)

; Draw the I0 and alpha numbers - these have no mag limits
newpar[0]=(pmax-pmin)[0]*randomu(seed)+pmin[0]
newpar[3]=(pmax-pmin)[3]*randomu(seed)+pmin[3]


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
      g=(pmax-pmin)[6:7]*randomu(seed,2)+pmin[6:7]
      newpar[6:7]=g
      gok=sqrt(total(g^2)) lt g_lim
   endif

   if not G1ok then begin
      G1=(pmax-pmin)[8:9]*randomu(seed,2)+pmin[8:9]
      newpar[8:9]=G1
      G1ok=sqrt(total(G1^2)) lt G1_lim
   endif

   if not G3ok then begin
      G3=(pmax-pmin)[8:9]*randomu(seed,2)+pmin[10:11]
      newpar[10:11]=G3
      G3ok=sqrt(total(G3^2)) lt G3_lim
   endif
   
   allok=ctrok and eok and gok and G1ok and G3ok

endwhile

return, newpar

end
