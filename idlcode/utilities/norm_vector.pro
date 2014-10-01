function norm_vector, V, ORIGIN=ORIGIN

dims=size(v,/dimensions)

if ((n_elements(dims) eq 1) or (dims[0] eq 1)) then return,sqrt(total(v^2))

vec=v

if (keyword_set(origin) and (n_elements(origin) eq 2)) then begin
   vec[*,0]-=origin[0]
   vec[*,1]-=origin[1]
endif

return,sqrt(total(vec^2,2))

end
