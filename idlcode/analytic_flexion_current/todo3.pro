pro todo3

  sel=1
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  if sel eq 1 then begin
  
    dir='data/a1689_hla/fourfilter/aim_fit/'
    fit=read_data(dir+'run4/run4_fit.dat')
    err=read_data(dir+'run4/run4_err.dat')
    fom=read_data(dir+'run4/run4_fom.dat')
    
    fill='expgrid500dev1000'
    ndev=1000
    ng=500L
    
    
    ok=where((fom[*,3] lt 1.5d) and (fom[*,11] eq 0) and $
      (fom[*,6] gt 1d-4) and (fom[*,6] lt 1d-2) and $
      (fit[*,6] gt 2d) and (fit[*,3]-2d*alog10(fit[*,6]) gt -1.5d),nok)
      
    x=dcomplex(fit[ok,1],fit[ok,2])
    ; psi1 = F/(1-kappa)
    psi1=dcomplex(fit[ok,11],fit[ok,12])*4d
    
    xmin=0d
    ymin=0d
    xmax=5800d
    ymax=5800d
    
    
    rf=[45d,60d,75d,90d]/0.05d
    lf=[3d,5d,7d,9d]
    head=dir+fill+'MKap'
    mid ='L'+['3','5','7','9']
    tail='R'+['45','60','75','90']
    
    
    window=mrdfits('data/a1689_hla/fourfilter/aim_fit/thesisv1/fieldwin.fits')
    
    ; Run through all the filters
    for i=0,n_elements(lf)-1 do for j=0,n_elements(rf)-1 do begin
    
      map=exp(-aperture_mass(x,psi1,lf[i],rf[j],nx=ng,ny=ng,/verbose,$
        xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax))
        
      ;Save the aperture mass map
      mwrfits,scale_image(map,5800d/ng)*window,head+'_'+mid[i]+'_'+tail[j]+'.fits'
      
      sig=map*0
      ave=map*0
      devs=dblarr(max([ndev,1]),ng,ng)
      
      if ndev gt 3 then begin
        print,'Deviating...'
        for k=0,ndev-1 do begin
          if ((k+1) mod 100) eq 0 then print,'On deviate #'+strcompress(string(k+1),/remove_all)
          psi1dev=psi1 + dcomplex(err[ok,9]*randomn(seed,nok),err[ok,10]*randomn(seed,nok))
          devs[k,*,*]=exp(-aperture_mass(x,psi1dev,lf[i],rf[j],nx=ng,ny=ng,$
            xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax))
        endfor
        
        for k=0,ng-1 do for m=0,ng-1 do begin
          sig[k,m]=stddev(devs[*,k,m])
          ave[k,m]=mean(devs[*,k,m])
        endfor
        
        zeroes=where((sig lt 1d-6) or (finite(sig) eq 0),nz)
        if nz gt 0 then sig[zeroes]=1d9
        
        
        mwrfits,scale_image(map/sig,5800d/ng)*window,head+'SNR_'+mid[i]+'_'+tail[j]+'.fits'
        mwrfits,scale_image(sig,5800d/ng)*window,head+'sig_'+mid[i]+'_'+tail[j]+'.fits'
        mwrfits,scale_image(ave,5800d/ng)*window,head+'ave_'+mid[i]+'_'+tail[j]+'.fits'
        
        print,'Done deviating!'
      endif
      
    endfor
    
    
    print,'Done!!!'
    
  endif
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  
  
end




