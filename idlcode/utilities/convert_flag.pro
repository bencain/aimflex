FUNCTION CONVERT_FLAG, FLAG, N_PAR, TO_2D=TO_2D, TO_1D=TO_1D

; This function converts an integer flag into a binary array or vice
; versa, depending on the TO_ND keywords.  If neither keyword is set,
; then the output is just the input.

if (keyword_set(to_1d) and keyword_set(to_2d)) then return,flag

if keyword_set(to_1d) then begin

   dims=size(flag,/dimensions)
   if n_elements(dims) eq 1 then useflag=transpose(flag) else useflag=flag
   dims=size(useflag,/dimensions)
   out=lonarr(dims[0])

   for i=0,dims[0]-1 do out[i]=long(total(2^double(where(useflag[i,*] eq 1))))

endif else if keyword_set(to_2d) then begin
   n_obj=n_elements(flag)
   out=lonarr(n_obj,n_par)
   for i=0,n_obj-1 do out[i,*]=long((reverse(binary(long(flag[i]))))[0:n_par-1])

   if n_obj eq 1 then out=transpose(out)
endif else out=flag

return,out

END
