function flexion_limit,lpars,G1lim=G1lim,G3lim=G3lim,glim=glim,count=count

; Assumes an Nx6 array of lensing parameters.  Will break if N
; isn't > 1.  These are limits on the magnitude of the fields

if not keyword_set(glim) then glim=1d3
if not keyword_set(G1lim) then G1lim=1d
if not keyword_set(G3lim) then G3lim=G1lim

g=sqrt(total(lpars[*,0:1]^2,2))
G1=sqrt(total(lpars[*,2:3]^2,2))
G3=sqrt(total(lpars[*,4:5]^2,2))

ok=where((g lt glim) and (G1 lt G1lim) and (G3 lt G3lim),count)

return,ok

end
