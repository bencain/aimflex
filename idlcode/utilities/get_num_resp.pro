FUNCTION GET_NUM_RESP, PROMPT, OK_RESP=OK_RESP, NTRIES=NTRIES, $
                       DEFAULT=DEFAULT

; This function prompts the user for a numerical response.  PROMPT is
; a string array which is printed to the screen.  The user is prompted
; NTRIES times for a valid response, after which a default value is
; returned. The valid responses are input in OK_RESP, but default to a
; binary 1/0 = yes/no answer.

if not keyword_set(ok_resp) then ok_resp=[0D,1D]
if not keyword_set(ntries) then ntries=10
if not keyword_set(default) then default=-1

ok_resp=long(ok_resp)
n_ok=n_elements(ok_resp)

print,'*******'
for i=0,n_elements(prompt)-1 do print,prompt[i]
print,'*******'

st='(Please enter '
for j=0,n_ok-2 do st+=$
   strcompress(string(ok_resp[j]),/remove_all)+', '
st+='or '+strcompress(string(ok_resp[n_ok-1]),/remove_all)+'): '

;print,st
resp=max(ok_resp)+5D

for i=0,ntries-1 do begin

   read,resp,prompt=st

   overlap=set_intersection(resp,ok_resp)

   if overlap eq -1 then begin
      badresp='Invalid response...'
      badresp+=strcompress('('+string(i+1)+'/'+string(ntries)+')',$
                           /remove_all)
      if i lt ntries-1 then badresp+=' Try again!' else $
         badresp+=' Defaulting to '+strcompress(string(default),/remove_all)
      print,badresp
   endif else return,resp[0]


endfor

return,default

end
