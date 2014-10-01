pro save_data_mesh, data_arr1, data_arr2, filename, comment=comment, $
                    delimiter=delimiter, tail=tail, lmesh=lmesh, rmesh=rmesh,$
                    format=format

; This procedure writes a 2-D array of data to a file given by
; filename.  The first line of the file will be a header line which
; has two elements: the sizes of the two dimensions of the array, in
; order. The second line will be a comment (blank if the comment
; keyword is not set).  The first two lines begin with a # sign.
; Following will be the elements of the 2-D array separated by the
; given delimiter and each line ending with the given tail.

if not keyword_set(delimiter) then delimiter=string(9B)
if not keyword_set(comment) then comment=''
if not keyword_set(tail) then tail=''
if not keyword_set(format) then format='(e11.4)'
if not keyword_set(lmesh) then lmesh='('
if not keyword_set(rmesh) then rmesh=')'


size1=size(data_arr1,/dimensions)
size2=size(data_arr2,/dimensions)
size=[min([size1[0],size2[0]]),min([size1[1],size2[1]])]

openw,wunit,filename,/get_lun,width=size[1]*17L

printf,wunit, print_prep(size,head='# ',format='(I7)',delimiter='  ')
printf,wunit,'# '+comment

a=dblarr(size[1])
b=dblarr(size[1])
for i=0,size[0]-1 do begin
   a[*]=data_arr1[i,*]
   b[*]=data_arr2[i,*]
   printf,wunit,print_prep_mesh(a,b,format=format,delimiter=delimiter,$
                                tail=tail,rmesh=rmesh,lmesh=lmesh)
endfor

close,wunit
free_lun,wunit

end



