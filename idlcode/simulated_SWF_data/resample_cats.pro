pro resample_cats, orig_fl_cat, orig_wl_cat, nresample


; Get the catalogs
;  print,'Reading WL'
  readcol,orig_wl_cat,wlra,wldec,wle1,wle2,wlwt,wlzed
;  print,'Reading FL'
  readcol,orig_fl_cat,flra,fldec,fle1,fle2,psi11,psi12,psi31,psi32,flwt,flzed
;  print,'Done reading'

  realformat='(E15.8)'

  nwl=n_elements(wlra)
  nfl=n_elements(flra)

; Resample the catalogs, with replacement
  for i=1,nresample do begin
     wl_sel=floor(double(nwl)*randomu(seed,nwl,/double),/L64)
     if nfl gt 0 then fl_sel=floor(double(nfl)*randomu(seed,nfl,/double),/L64)

; Sort out the selections in increasing order, find unique selections
     wl_sel=wl_sel[sort(wl_sel)]
     if nfl gt 0 then fl_sel=fl_sel[sort(fl_sel)]

     wl_out=wl_sel[uniq(wl_sel)]
     if nfl gt 0 then fl_out=fl_sel[uniq(fl_sel)]

     nwl_out=n_elements(wl_out)
     if nfl gt 0 then nfl_out=n_elements(fl_out) else nfl_out=0L

; Deal with the repeats by increasing the weights
     wl_out_wt=wlwt[wl_out]
     if nfl gt 0 then fl_out_wt=flwt[fl_out]

     for j=0,nwl_out-1 do wl_out_wt[j]*=double(n_elements(where(wl_sel eq wl_out[j])))
     for j=0,nfl_out-1 do fl_out_wt[j]*=double(n_elements(where(fl_sel eq fl_out[j])))

     wlfile=file_dirname(orig_wl_cat,/mark_directory)+'wl_tmp'+string(i,format='(I03)')+'.cat'
     flfile=file_dirname(orig_fl_cat,/mark_directory)+'fl_tmp'+string(i,format='(I03)')+'.cat'
     openw,wlunit,wlfile,/get_lun,width=200
     openw,flunit,flfile,/get_lun,width=200
     for j=0,nwl_out-1 do printf,wlunit,$
                                 string(wlra[wl_out[j]],format=realformat)+' '+$
                                 string(wldec[wl_out[j]],format=realformat)+' '+$
                                 string(wle1[wl_out[j]],format=realformat)+' '+$
                                 string(wle2[wl_out[j]],format=realformat)+' '+$
                                 string(wl_out_wt[j],format=realformat)+' '+$
                                 string(wlzed[wl_out[j]],format=realformat)
     for j=0,nfl_out-1 do printf,flunit,$
                                 string(flra[fl_out[j]],format=realformat)+' '+$
                                 string(fldec[fl_out[j]],format=realformat)+' '+$
                                 string(fle1[fl_out[j]],format=realformat)+' '+$
                                 string(fle2[fl_out[j]],format=realformat)+' '+$
                                 string(psi11[fl_out[j]],format=realformat)+' '+$
                                 string(psi12[fl_out[j]],format=realformat)+' '+$
                                 string(psi31[fl_out[j]],format=realformat)+' '+$
                                 string(psi32[fl_out[j]],format=realformat)+' '+$
                                 string(fl_out_wt[j],format=realformat)+' '+$
                                 string(flzed[fl_out[j]],format=realformat)
     close,wlunit,flunit
     free_lun,wlunit,flunit

  endfor



end
