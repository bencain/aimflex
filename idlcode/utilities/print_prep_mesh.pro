function print_prep_mesh,array1,array2, head=head, delimiter=delimiter, $
                         tail=tail, lmesh=lmesh, rmesh=rmesh, format=format
  
; This function takes two 1-D array and returns a string of the
; following format:

; HEAD A1_1 LMESH A2_1 RMESH DELIMITER A1_2 LMESH A2_2 RMESH DELIMITER ...
; ... DELIMITER A1_N LMESH A2_N RMESH TAIL

;where N is the minimum of the number of elements in array1 and
;array2.


tab=string(9B)

if not keyword_set(delimiter) then delimiter=tab
if not keyword_set(format) then format='(e11.4)'
if not keyword_set(lmesh) then lmesh='('
if not keyword_set(rmesh) then rmesh=')'

n=min([n_elements(array1),n_elements(array2)])

if keyword_set(head) then begin
   out=head+strcompress(string(array1[0],format=format),/remove_all)+lmesh+$
       strcompress(string(array2[0],format=format),/remove_all)+rmesh+$
       delimiter
endif else begin
   out=strcompress(string(array1[0],format=format),/remove_all)+lmesh+$
       strcompress(string(array2[0],format=format),/remove_all)+rmesh+$
       delimiter
endelse

for i=1,n-2 do begin
   out+=strcompress(string(array1[i],format=format),/remove_all)+lmesh+$
        strcompress(string(array2[i],format=format),/remove_all)+rmesh+$
        delimiter
endfor

out+=strcompress(string(array1[n-1],format=format),/remove_all)+lmesh+$
     strcompress(string(array2[n-1],format=format),/remove_all)+rmesh


if keyword_set(tail) then out+=tail

return,out
end

