function read_data, filename, comment=comment, size=size, quiet=quiet

; This function reads a 2-D array from the file indicated by filename
; and returns it.  The file format is that specified for save_data.
; The comment is placed into the comment keyword and the dimensions of
; the data array is put into size.

openr,runit,filename,/get_lun

;Get the data array size
size_string=''
readf,runit,size_string
ss=strparse(size_string)
size=long(ss[1:2])

;Get the comment and remove the leading '#" character
comment=string(1)
readf,runit,comment

comment=byte(comment)
nchar=n_elements(comment)
if nchar gt 3 then comment=string(comment[2:nchar-1]) else comment=''

;Get the data
data=dblarr(size)
temp=dblarr(size[1])

for i=0,size[0]-1 do begin
   readf,runit,temp
   data[i,*]=double(temp)
endfor

if not keyword_set(quiet) then begin
   print,filename+' Comment:'
   print,comment
endif

close,runit,/all

return,data

end
