pro zedshuffle, orig_fl_cat, orig_wl_cat, dzphot, ncat


; Get the catalogs
;  print,'Reading WL'
  readcol,orig_wl_cat,wlra,wldec,wle1,wle2,wlwt,wlzed
;  print,'Reading FL'
  readcol,orig_fl_cat,flra,fldec,fle1,fle2,psi11,psi12,psi31,psi32,flwt,flzed
;  print,'Done reading'
  
  realformat='(E15.8)'

  nwl=n_elements(wlra)
  nfl=n_elements(flra)
  
; The error is usually characterized in terms of scatter in dz/(1+z).
; We need to calculate dz
  wl_err=(1d + wlzed)*dzphot
  fl_err=(1d + flzed)*dzphot

  for i=1,ncat do begin

; Calculate the new redshifts
     wlzed_new=wlzed + wl_err*randomn(seed,nwl,/double)
     flzed_new=flzed + fl_err*randomn(seed,nfl,/double)

; Name and open files
     wlfile=file_dirname(orig_wl_cat,/mark_directory)+'wl_tmp'+string(i,format='(I03)')+'.cat'
     flfile=file_dirname(orig_fl_cat,/mark_directory)+'fl_tmp'+string(i,format='(I03)')+'.cat'

     openw,wlunit,wlfile,/get_lun,width=200
     openw,flunit,flfile,/get_lun,width=200



     for j=0,nwl-1 do printf,wlunit,$
                             string(wlra[j],format=realformat)+' '+$
                             string(wldec[j],format=realformat)+' '+$
                             string(wle1[j],format=realformat)+' '+$
                             string(wle2[j],format=realformat)+' '+$
                             string(wlwt[j],format=realformat)+' '+$
                             string(wlzed_new[j],format=realformat)
     for j=0,nfl-1 do printf,flunit,$
                             string(flra[j],format=realformat)+' '+$
                             string(fldec[j],format=realformat)+' '+$
                             string(fle1[j],format=realformat)+' '+$
                             string(fle2[j],format=realformat)+' '+$
                             string(psi11[j],format=realformat)+' '+$
                             string(psi12[j],format=realformat)+' '+$
                             string(psi31[j],format=realformat)+' '+$
                             string(psi32[j],format=realformat)+' '+$
                             string(flwt[j],format=realformat)+' '+$
                             string(flzed_new[j],format=realformat)
     close,wlunit,flunit
     free_lun,wlunit,flunit
     
  endfor



end
