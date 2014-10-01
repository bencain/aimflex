function density_map, pos, nlimit=nlimit

if not keyword_set(nlimit) then nlimit=6

dims=size(pos,/dimensions)

if n_elements(dims) ne 2 then return,0d
if dims[1] lt 2 then return,dims[*,0]*0d
if dims[0] lt nlimit then return,dims[*,0]*0d

sigma=dblarr(dims[0])

for i=0,dims[0]-1 do begin

   rs=sqrt((pos[*,0]-pos[i,0])^2+(pos[*,1]-pos[i,1])^2)

   use=(sort(r))[0:nlimit]

   r0=max(rs[use])

   sigma[i]=double(nlimit+1)/(!dpi*r0^2)

endfor

return,density

end
