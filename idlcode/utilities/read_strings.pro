function read_strings, filename, comment=comment, size=size, quiet=quiet

openr,runit,filename,/get_lun

;Get the data array size
size_string=''
readf,runit,size_string
ss=strparse(size_string)
size=long(ss[1])

;Get the comment and remove the leading '#" character
comment=string(1)
readf,runit,comment

comment=byte(comment)
nchar=n_elements(comment)
if nchar gt 2 then comment=string(comment[2:nchar-1]) else comment=''

;Get the data
strings=string(dblarr(size))
temp=' '

for i=0,size[0]-1 do begin
   readf,runit,temp
   strings[i]=temp
endfor

if not keyword_set(quiet) then begin
   print,filename+' Comment:'
   print,comment
endif

return,strings

end
