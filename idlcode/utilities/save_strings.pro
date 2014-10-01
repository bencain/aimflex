pro save_strings, string_arr, filename, comment=comment

; This procedure writes a 2-D array of data to a file given by
; filename.  The first line of the file will be a header line which
; has two elements: the sizes of the two dimensions of the array, in
; order. The second line will be a comment (blank if the comment
; keyword is not set).  The first two lines begin with a # sign.
; Following will be the elements of the 2-D array separated by the
; given delimiter and each line ending with the given tail.

if not keyword_set(comment) then comment=''

size=n_elements(string_arr)

openw,wunit,filename,/get_lun

printf,wunit, print_prep(size,head='# ',format='(I7)',delimiter='  ')
printf,wunit,'# '+comment

if finite(string_arr[0]) then begin
   for i=0,size-1 do begin
      printf,wunit,string_arr[i]
   endfor
endif else begin
   printf,wunit,'*** NAN ***'
endelse

close,wunit
free_lun,wunit

end



