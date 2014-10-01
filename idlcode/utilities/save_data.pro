pro save_data, data_arr, filename, comment=comment, delimiter=delimiter, $
               tail=tail, format=format, noheader=noheader, nocomment=nocomment

; This procedure writes a 2-D array of data to a file given by
; filename.  The first line of the file will be a header line which
; has two elements: the sizes of the two dimensions of the array, in
; order. The second line will be a comment (blank if the comment
; keyword is not set).  The first two lines begin with a # sign.
; Following will be the elements of the 2-D array separated by the
; given delimiter and each line ending with the given tail.  If the
; array is 1-D, then it will be transposed to a 1 x N array, and if
; the array is k-D with k > 2 only the first two dimensions will be
; considered. 

if not keyword_set(delimiter) then delimiter=string(9B)
if not keyword_set(comment) then comment=''
if not keyword_set(tail) then tail=' '
if not keyword_set(format) then format='(e18.10)'

ndim=n_elements(size(data_arr,/dimensions))

if ndim lt 2 then begin
   arr=transpose(data_arr)
endif else arr=data_arr
size=size(arr,/dimensions)

openw,wunit,filename,/get_lun,width=max([size[1]*25L,25L])

if not keyword_set(noheader) then begin
   printf,wunit, print_prep(size,head='# ',format='(I7)',delimiter='  ',$
                            tail=' ')
endif
if not keyword_set(nocomment) then begin
   ncom=n_elements(comment)
   for k=0,ncom-1 do printf,wunit,'# '+comment[k]
endif

if finite(data_arr[0]) then begin

   a=dblarr(size[1])

   for i=0,size[0]-1 do begin
      a[*]=arr[i,*]
      printf,wunit,print_prep(a,format=format,delimiter=delimiter,tail=tail)
   endfor
endif else begin
   printf,wunit,'*** NAN ***'
endelse

close,wunit
free_lun,wunit

end



