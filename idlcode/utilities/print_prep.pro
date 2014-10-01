function print_prep,array, head=head, delimiter=delimiter, tail=tail, $
                    format=format

; This function takes a 1-D array and returns a string of each of the
; elements with the given delimiter between.  The delimiter is taken
; to be a tab if not set.  If the head keyword is set, then the output
; string begins with head, and if the tail keyword is set then the
; string ends with tail.

tab=string(9B)

if not keyword_set(delimiter) then delimiter=tab
if not keyword_set(format) then format='(e11.4)'

n=n_elements(array)

if keyword_set(head) then begin
   out=head+strcompress(string(array[0],format=format),/remove_all)+$
       delimiter
endif else begin
   out=strcompress(string(array[0],format=format),/remove_all)+delimiter
endelse

for i=1,n-2 do out+=strcompress(string(array[i],format=format),/remove_all)+$
   delimiter

out+=strcompress(string(array[n-1],format=format),/remove_all)


if keyword_set(tail) then out+=tail

return,out
end

