pro bootstrap_medianstats, indir, nfile, label, field, outdir

  files=indir+'/run.'+label+'nofl_'+string(lindgen(nfile),format='(I03)')+'_'+field+'_rs_inf.fits'

  img0=mrdfits(files[0],0,hdr)
  npix=size(img0,/dimensions)

  dims=[nfile,npix[0],npix[1]]

  imgs=dblarr(dims)

  imgs[0,*,*]=img0

  for i=1,nfile-1 do begin
     imgs[i,*,*]=mrdfits(files[i])
  end
  
  img16=dblarr(npix)
  img50=dblarr(npix)
  img84=dblarr(npix)
  
  for i=0,npix[0]-1 do for j=0,npix[1]-1 do begin
     
     pcts=percentiles(imgs[*,i,j],value=[0d,0.16d,0.5d,0.84d,1d])
     img16[i,j]=pcts[1]
     img50[i,j]=pcts[2]
     img84[i,j]=pcts[3]


  end

  mwrfits,img16,outdir+'/'+label+'_'+field+'_16pct.fits',hdr
  mwrfits,img50,outdir+'/'+label+'_'+field+'_50pct.fits',hdr
  mwrfits,img84,outdir+'/'+label+'_'+field+'_84pct.fits',hdr
  
  
end

