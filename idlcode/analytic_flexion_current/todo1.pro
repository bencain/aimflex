pro todo1                                                     

sel=4

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if sel eq 1 then begin

dir='data/a1689_hla/fourfilter/aim_fit/'
fit=read_data(dir+'run4/run4_fit.dat')
err=read_data(dir+'run4/run4_err.dat')
fom=read_data(dir+'run4/run4_fom.dat')

ok=where((fom[*,3] lt 1.5d) and (fom[*,11] eq 0) and $
         (fom[*,6] gt 1d-4) and (fom[*,6] lt 1d-2) and $
         (fit[*,6] gt 2d) and (fit[*,3]-2d*alog10(fit[*,6]) gt -1.5d),nok)


fill='abs300_'
x=dcomplex(fit[ok,1],fit[ok,2])
psi1=dcomplex(fit[ok,11],fit[ok,12])

xmin=0d
ymin=0d
xmax=5800d
ymax=5800d


rf=[45d,60d,75d,90d]/0.05d
lf=[3d,5d,7d];,9d]
head=dir+fill+'MKap'
mid ='L'+['3','5','7'];,'9']
tail='R'+$
     ['45','60','75','90']

ndev=100
ng=300L


window=mrdfits('data/a1689_hla/fourfilter/aim_fit/thesisv1/fieldwin.fits')

for i=0,n_elements(lf)-1 do for j=0,n_elements(rf)-1 do begin

   map=aperture_mass(x,psi1,lf[i],rf[j],nx=ng,ny=ng,/verbose,$
                     xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)
;   map=abs(aperture_mass(x,psi1,lf[i],rf[j],nx=ng,ny=ng,/verbose,$
;                         xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax))
   
   
   sig=map*0
   ave=map*0
   devs=dblarr(max([ndev,1]),ng,ng)

   print,'Deviating...'
   for k=0,ndev-1 do begin
      psi1dev=psi1 + dcomplex(err[ok,9]*randomn(seed,nok),err[ok,10]*randomn(seed,nok))
      devs[k,*,*]=aperture_mass(x,psi1dev,lf[i],rf[j],nx=ng,ny=ng,$
                                xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)
;      devs[k,*,*]=abs(aperture_mass(x,psi1dev,lf[i],rf[j],nx=ng,ny=ng,$
;                                    xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax))

   endfor

   for k=0,ng-1 do for m=0,ng-1 do begin
      sig[k,m]=stddev(devs[*,k,m])
   endfor
   
   zeroes=where((sig lt 1d-6) or (finite(sig) eq 0),nz)
   if nz gt 0 then sig[zeroes]=1d9


;   mwrfits,scale_image(map/sig,5800d/ng)*window,head+'SNR_'+mid[i]+'_'+tail[j]+'.fits'
   mwrfits,scale_image(abs(map)/sig,5800d/ng)*window,head+'SNR_'+mid[i]+'_'+tail[j]+'.fits'

   print,'Done deviating!'
endfor

print,'Done!!!'

endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


if sel eq 3 then begin

snrs=['gopsiSNR_L3_R60.fits',$
      'gopsiSNR_L3_R75.fits',$
      'gopsiSNR_L5_R60.fits',$
      'gopsiSNR_L5_R75.fits',$
      'gopsiSNR_L7_R60.fits',$
      'gopsiSNR_L7_R75.fits']

dir='data/a1689_hla/fourfilter/aim_fit/'

field=mrdfits(dir+'fieldwin.fits')

print,max(field),min(field)

for i=0,n_elements(snrs)-1 do begin
   im=mrdfits(dir+snrs[i])
   win=im*0d + 1d
   ch1=3d*aim_mk_image([3.9d,193d - 250d,312d - 250d,37d,0d,0.1d,dblarr(6)],0d,win,/nopsf)
   ch2=3d*aim_mk_image([3.9d,263d - 250d,211d - 250d,37d,0d,0.1d,dblarr(6)],0d,win,/nopsf)
   
   im=scale_image(im+ch1+ch2,11.6)
   im*=field
   mwrfits,im,dir+'ok'+snrs[i]
endfor

endif

if sel eq 4 then begin
   dir='data/a1689_hla/fourfilter/aim_fit/'
   radec=read_data(dir+'use_radec.dat')
   fit=read_data(dir+'use_fit.dat')
   err=read_data(dir+'use_err.dat')

   ok=where((err[*,9]/0.05d gt 1d-4) and (err[*,10]/0.05d gt 1d-4) and $
            (err[*,11]/0.05d gt 1d-4) and (err[*,12]/0.05d gt 1d-4))
   radec=radec[ok,*]
   fit=fit[ok,*]
   err=err[ok,*]
   ord=(sort(sqrt(err[*,9]^2+err[*,10]^2)/sqrt(fit[*,11]^2+fit[*,12]^2)))[0:99]

   radec=radec[ord,*]
   fit=fit[ord,*]
   err=err[ord,*]

   snr=sqrt(fit[*,11]^2+fit[*,12]^2)/sqrt(err[*,9]^2+err[*,10]^2)


   rah=radec[*,1]/15d
   hrs=floor(rah)
   min=floor((rah-hrs)*60d)
   sec=floor(((rah-hrs)*60d - min)*6000d)/100d

   absdec=abs(radec[*,2])
   signdec=radec[*,2]/absdec
   deg=floor(absdec)
   am=floor((absdec-deg)*60d)
   as=floor(((absdec-deg)*60d - am)*600d)/10d
   

   openw,shapeu1,dir+'tab_shape1.txt',/get_lun
   for i=0,99 do $
      printf,shapeu1,$
             i+1,hrs[i],min[i],sec[i],signdec[i]*deg[i],am[i],as[i],$ ; N, Positions
             fit[i,3],max([err[i,1],1d-4]),$                          ; logS0
             fit[i,4]*0.05d,max([err[i,2],2d-3])*0.05d,$              ; Xc
             fit[i,5]*0.05d,max([err[i,3],2d-3])*0.05d,$              ; Yc
             format='(I3," & ",I2,":",I2,":",F05.2," & ",I3,":",I2,":",F04.1," & ",'+$ ; N, RA, Dec
             'F7.3," & ",F7.4," & ",F7.3," & ",F7.4," & ",F7.3," & ",F7.4," \\")'      ; logS0, Xc, Yc
   close,shapeu1
   openw,shapeu2,dir+'tab_shape2.txt',/get_lun
   for i=0,99 do $
      printf,shapeu2,$
             i+1,fit[i,6]*0.05d,max([err[i,4],2d-3])*0.05d,$ ; N, alpha
             fit[i,7],max([err[i,5],1d-4]),$     ; E+
             fit[i,8],max([err[i,6],1d-4]),$     ; Ex
             format='(I3," & ",F7.4," & ",F7.4," & ",F7.3," & ",F7.4," & ",F7.3," & ",F7.4," \\")'
   close,shapeu2


   openw,flexu,dir+'tab_flex.txt',/get_lun
   for i=0,99 do $
      printf,flexu,$
             i+1,$
             fit[i,11]/0.05d,err[i,9]/0.05d,fit[i,12]/0.05d,err[i,10]/0.05d,$  ; Psi11, Psi12
             fit[i,13]/0.05d,err[i,11]/0.05d,fit[i,14]/0.05d,err[i,12]/0.05d,$ ; Psi31, Psi32
;             format='(I3," & ",F7.4," & ",F7.4," & ",F7.4," & ",F7.4," & ",'+$ ; N, Psi11, Psi12
;             'F7.4," & ",F7.4," & ",F7.4," & ",F7.4," \\")'                    ; Psi31, Psi32
             snr[i],$                                                          ; SNR
             format='(I3," & ",F7.4," & ",F7.4," & ",F7.4," & ",F7.4," & ",'+$ ; N, Psi11, Psi12
             'F7.4," & ",F7.4," & ",F7.4," & ",F7.4," & ",F5.1," \\")'         ; Psi31, Psi32, SNR
   close,flexu

   corr=read_data(dir+'use_corr.dat')

   print,ord[5]
   print,convert_sym_matrix(corr[ord[5],1:*])
   print,err[5,1:12]*[1d,0.05d,0.05d,0.05d,1d,1d,1d,1d,20d,20d,20d,20d]

endif

if sel eq 5 then begin
   i=-1 + 5
   dir='data/20101102/'
   fitfiles=dir+'run0'+['1','2','3','4','5','6']+'_fit.dat'
   fomfiles=dir+'run0'+['1','2','3','4','5','6']+'_fom.dat'
   trufiles=dir+'run0'+['1','2','3','4','5','6']+'_true.dat'
   stafiles=dir+'run0'+['1','2','3','4','5','6']+'_start.dat'
   errfiles=dir+'run0'+['1','2','3','4','5','6']+'_err.dat'
   tags='plot'+['1','2','3','4','5','6']

   for i=2,5 do $
      aim_mcmc_plots, fitfiles[i],fomfiles[i],trufiles[i],errfiles[i],stafiles[i],tag=tags[i]

endif


end




