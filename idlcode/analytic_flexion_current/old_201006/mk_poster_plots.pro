pro mk_poster_plots,fitfile,chisqcut=chisqcut,G1cut=G1cut,G3cut=G3cut

fit=read_data(fitfile)

chisq=fit[*,17]

g=fit[*,9:10]
G1=fit[*,11:12]
G3=fit[*,13:14]

ok=lindgen(n_elements(chisq))

; Make a cut on chisq/dof first
if not keyword_set(chisqcut) then chisqcut=1d

chisq_ok=where(chisq lt chisqcut)
ok=set_intersection(ok,chisq_ok)

; Now make cuts based on the lensing parameters
if not keyword_set(G1cut) then G1cut=3d
if not keyword_set(G3cut) then G3cut=3d

G1_ok=where(sqrt(total(G1^2,2)) lt G1cut)
G3_ok=where(sqrt(total(G3^2,2)) lt G3cut)

Gn_ok=set_intersection(G1_ok,G3_ok)

ok=set_intersection(ok,Gn_ok)

;plot_field,fit[ok,*],xy_cols=[1,2],field_cols=[11,12]
radial_flexvar,fit[ok,*],[2500d,2500d],nbin=8,nwedge=1,nobj=nobj;,$
;               pstag='test',/to_eps,ndensity=ndensity



end
