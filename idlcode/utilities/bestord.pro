pro bestord,nlin,blaz,lam,ord,ang,eff,ord1=ord1,nord=nord,$
            minang=minang,maxang=maxang,help=help
;procedure to calculate the best order for a grating
;at a given wavelength and to give the angle
;INPUTS:
;  nlin - number of lines per mm of the grating
;  blaz - the blaze angle of the grating [degrees]
;  lam - wavelengths to consider [microns]
;OUTPUTS:
;  ord - order that puts LAM closest to peak of blaze function
;        NOTE: resolution desired may cause user to choose different
;              order.  The higher the order, the steeper the grating
;              angle, the higher the resolution.
;  ang - the grating angle for given wavelength [degrees]
;  eff - the efficiency of the grating at a given grating angle
;          calculated assuming each order's efficiency curve is
;          essentially a Gaussian with FWHM equal to the order spacing
;KEYWORDS:
;  ord1=ord1 - lowest order expected [default is 1]
;  nord=nord - number of orders to consider [default is 30]
;  minang=minang - minimum angle [degrees] allowed for grating [default is 2]
;  maxang=maxang - maximum angle [degrees] allowed for grating [default is 75]
;
;mjr Jun-1997: create

if keyword_set(help) then begin
  print,'BESTORD,nlin,blaz,lam,ord[,ang,eff][,ord1=ord1 (def=1)],$'
  print,'        [,nord=nord (def=30)][,minang=minang (def=2)],$'
  print,'        [,maxang=maxang (def=75)][,/help]'
  print,'procedure to calculate the best order for a grating'
  print,'at a given wavelength and to give the angle'
  print,'INPUTS:'
  print,'  nlin - number of lines per mm of the grating'
  print,'  blaz - the blaze angle of the grating [degrees]'
  print,'  lam - wavelengths to consider [microns]'
  print,'OUTPUTS:'
  print,'  ord - order that puts LAM closest to peak of blaze function'
  print,'        NOTE: resolution desired may cause user to choose different'
  print,'              order.  The higher the order, the steeper the grating'
  print,'              angle, the higher the resolution.'
  print,'  ang - the grating angle for given wavelength [degrees]'
  print,'  eff - the efficiency of the grating at a given grating angle'
  print,'          calculated assuming each order''s efficiency curve is'
  print,'          essentially a Gaussian with FWHM equal to the order spacing'
  print,'KEYWORDS:'
  print,'  ord1=ord1 - lowest order expected [default is 1]'
  print,'  nord=nord - number of orders to consider [default is 30]'
  print,'  minang=minang - minimum angle [degrees] allowed for '
  print,'                  grating [default is 2]'
  print,'  maxang=maxang - maximum angle [degrees] allowed for '
  print,'                  grating [default is 75]'
  retall
endif

if n_params() lt 4 then begin
  print,'BESTORD,nlin,blaz,lam,ord[,ang,eff][,ord1=ord1 (def=1)],$'
  print,'        [,nord=nord (def=30)][,minang=minang (def=2)],$'
  print,'        [,maxang=maxang (def=75)][,/help]'
  retall
endif

  if max(lam) ge 2*1d3/nlin then begin
    print,'WARNING: groove spacing is shorter than longest wavelength'
    print,'groove spacing [micron] = '+strtrim(1d3/nlin,2)
    retall
  endif

  if not(keyword_set(ord1)) then ord1=1
  if not(keyword_set(nord)) then nord=30
  if not(keyword_set(minang)) then minang=2
  if not(keyword_set(maxang)) then maxang=75

  nlam=n_elements(lam)		;number of wavelength points
  allord=indgen(nord)+ord1	;ALL ORDers considered
  r_d=18d1/!pi			;convert from Radians to Degrees
  tblaz=blaz/r_d		;BLAZe angle [radians]
  blam=2*1d3/nlin*sin(tblaz)/allord	;calculate blaze wavelengths
  ord=intarr(nlam)
  ang=dblarr(nlam)
  eff=dblarr(nlam)
  for i=0,nlam-1 do begin
    dum=min((blam-lam(i))^2,m)
    if lam(i)*allord(m)/2*nlin/1d3 gt sin(maxang/r_d) then m=m-1>1
    if lam(i)*allord(m)/2*nlin/1d3 lt sin(minang/r_d) then m=m+1
    ord1=allord(m)
    ord2=allord(m-1>0)
    eff1=exp(-.5*(blam(m)-lam(i))^2/(blam(m)/ord1/2.35)^2)
    eff2=exp(-.5*(blam(m-1>0)-lam(i))^2/(blam(m-1>0)/ord2/2.35)^2)
    if eff1 ge eff2 then begin
      ord(i)=ord1
      eff(i)=eff1
    endif else begin
      ord(i)=ord2
      eff(i)=eff2
    endelse
  endfor
  ang=asin(lam*ord/2*nlin/1d3)*r_d
end

