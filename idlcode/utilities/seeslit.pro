pro seeslit,slitwid,slitlen,sliteff
;procedure to estimate the slit efficiency when the
;spot at the slit is governed by seeing (Gaussian)
;INPUTS:
;  slitwid - width of slit in FWHM of the seeing disk
;  slitlen - length of slit to consider in FWHM of seeing disk
;OUTPUTS:
;  sliteff
;KEYWORDS:
;
;mjr Jun-1997: create

if n_params() lt 3 then begin
  print,'SEESLIT,slitwid,slitlen,sliteff'
  retall
endif

  common seemod,seex,seeing

  sz=size(seex)
  if sz(1) ne 1000 then begin
    dist_circle,seex,[1000,1000],499.5,499.5
    seex=seex/100.
    seeing=exp(-(seex^2/2./alog(2.)))
  endif

  dum=where(seex(*,500) le slitwid)
  dumx=minmax(dum)
  dum=where(seex(*,500) le slitlen)
  dumy=minmax(dum)
  sliteff=total(seeing(dumx(0):dumx(1),dumy(0):dumy(1)))/total(seeing)

end
