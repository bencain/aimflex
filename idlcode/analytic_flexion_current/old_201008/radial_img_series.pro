PRO RADIAL_IMG_SERIES,EPARS=EPARS,NGAL=NGAL,THETA_E=THETA_E,TAG=TAG


if not keyword_set(epars) then epars=[15d,0d,0d]
if not keyword_set(ngal) then ngal=20
if not keyword_set(theta_e) then theta_e=200d
if not keyword_set(tag) then tag='radial_im_series_' else tag+='_'


rgal=theta_e*(0.55001d +0.45d*dindgen(ngal)/(ngal-1))

pars=dblarr(12)
pars[0]=0.3d
pars[3:5]=epars

bkg=1000d
win=mk_window(5d*1.7d*epars[0]+1d)
psf=mk_gaussian_psf(0.25d)

for i=0,ngal-1 do begin

   fitsname=tag+'im'

   if i lt 9 then fitsname+='0'
   if i lt 99 then fitsname+='0'
   fitsname+=strcompress(string(i+1),/remove_all)+'.fits'

   lpars=sis_lensmodel([rgal[i],0],theta_e=theta_e)
   pars[6:11]=lpars

   print,[rgal[i],lpars[[0,2,4]]]

   im=mk_image(pars,bkg,win,psf)

   mwrfits,im,fitsname

endfor


end
