pro parameter_status,pars,pnames=pnames,wnum,label=label,delete=delete

if ((n_elements(wnum) eq 0) or  (!d.window eq -1)) then begin
   window,xsize=350,ysize=500,/free,title= '  '
   wnum=!d.window
endif else begin
   wset,wnum
endelse

npar=n_elements(pars)
if ((not keyword_set(pnames)) or (n_elements(pnames) ne npar)) then begin
   pnames=strarr(npar)
   for i=0,npar-1 do begin
      if i lt 9 then fill=' ' else fill=''
      pnames[i]='Par '+fill+strcompress(string(i+1),/remove_all)
   endfor
endif

if keyword_set(delete) then wdelete,wnum else begin

   if not keyword_set(label) then label='*** PARAMETER STATUS ***'

   plot,[0,1],[0,0],ystyle=4,xstyle=4,ytitle=label,$
        charsize=2d,xmargin=[20,5],ymargin=[0,0]
   ; Now output the text
   ys=(475d - (450d/double(npar-1))*dindgen(npar))/500d

   for i=0,npar-1 do begin
      xyouts,60d/350d,ys[i],pnames[i]+' = '+$
             strcompress(string(pars[i]),/remove_all),$
             /norm,charsize=2d,charthick=2d
   endfor
   xyouts,30d/350d,0.25d,label,/norm,charsize=2d,charthick=2d,orientation=90
endelse


end
