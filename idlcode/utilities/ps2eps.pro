pro PS2EPS, FILE, DELETE=DELETE, DELOK=DELOK

; This procedure rotates and fixes postscript files by using a series
; of commands sent to the command line OS via SPAWN

delok=1
split=strsplit(file,'.',/extract)
n=n_elements(split)

if n gt 1 then begin
   if (split[n-1] eq 'ps') then begin
      epsfile=''
      for i=0,n-2 do epsfile+=split[i]+'.'
      epsfile+='eps'
      if (file_test(epsfile) eq 0) then begin
         spawn,'/usr/texbin/ps2eps -B -q '+file
         if keyword_set(delete) then spawn,'rm -f '+file
      endif else print,'PS2EPS.PRO Error: EPS file already exists!'
   endif else print,'PS2EPS.PRO Error: Bad file name!'
endif else print,'PS2EPS.PRO Error: Bad file name!'

delok=file_test(file) and keyword_set(delete)

end



