function determ, M

s=size(M)
kill=0

if s[0] ne 2 then begin
   print,'DETERM.PRO -- Only works for 2-D matrices.'
   kill=1
endif

if s[1] ne s[2] then begin
   print,'DETERM.PRO -- Non-square matrix alert!'
   kill=1
endif

type=s[3]

if kill eq 1 then return,(make_array(1,type=type,value=0))[0]

n=s[1]

neg=(make_array(1,type=type,value=-1))[0]

if n eq 2 then begin

   det=M[0,0]*M[1,1]-M[0,1]*M[1,0]

endif else begin

   det=(make_array(1,type=type,value=0))[0]

   for i=0L,n-1 do begin

      cols=set_difference(lindgen(n),i)
      rows=lindgen(n-1)+1L

      subM=(M[cols,*])[*,rows]

      det+=M[i,0]*determ(subm)*neg^i

   endfor
endelse

return,det

end
