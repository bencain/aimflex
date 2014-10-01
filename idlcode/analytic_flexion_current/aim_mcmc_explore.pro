pro aim_mcmc_explore, n_mcmc, img_noise, tag, $
                     fixshear=fixshear, fixcenter=fixcenter,$
                     psffwhm=psffwhm,shearnoise=shearnoise, $
                     detJaclim=detJaclim, maxiter=maxiter, $
                     cutrad_scale=cutrad_scale, $
                     starttrueflex=starttrueflex,$
                     sis_lens=sis_lens, nfw_lens=nfw_lens, nfw_scale=nfw_scale,$
                     nsave=nsave, const_error=const_error,$
                     no_recenter_dataimg=no_recenter_dataimg,$
                     sersic_true=sersic_true, sersic_fit=sersic_fit,$
                     pseudogauss_true=pseudogauss_true, pseudogauss_fit=pseudogauss_fit,$
                     moffat_true=moffat_true, moffat_fit=moffat_fit,$
                     fixk1=fixk1, fixk2=fixk2, nopsi3=nopsi3

; This procedure makes a set of random input parameter values, creates
; a data image with the corresponding Gaussian noise properties and
; fits for the model parameters.

starttime=systime()

if img_noise le 0d then img_noise=0d
if not keyword_set(detJaclim) then detJaclim=0.2d
if keyword_set(shearnoise) then fixshear=1 else shearnoise=0d
if not keyword_set(psffwhm) then psffwhm=1.8d
if (keyword_set(sersic_true) or keyword_set(moffat_true)) then npar_true=13 else $
   if keyword_set(pseudogauss_true) then npar_true=14 else npar_true=12
if (keyword_set(sersic_fit) or keyword_set(moffat_fit)) then npar_fit=13 else $
   if keyword_set(pseudogauss_fit) then npar_fit=14 else npar_fit=12
if not keyword_set(maxiter) then maxiter=1000L

if not keyword_set(nsave) then nsave=-1

startps=dblarr(n_mcmc,npar_fit+1) ; An index number plus starting parameters
trueps=dblarr(n_mcmc,npar_true+1) ; An index number plus true parameters
fitps=dblarr(n_mcmc,npar_fit+1)   ; An index number plus fit parameters
errps=dblarr(n_mcmc,npar_fit+1)   ; An index number plus parameter errors
foms=dblarr(n_mcmc,12)            ; An index number plus figures of merit:
                                  ;    chisq, DoF, chisq/DoF, D1, D3,
                                  ;    psi1_err, psi3_err, 
                                  ;    pegflag, niter, failflag, badflag
nfit=npar_fit
if keyword_set(fixshear) then nfit-=2
if keyword_set(fixcenter) then nfit-=2
corrs=dblarr(n_mcmc,(nfit^2-nfit)/2+1) ; A place to save the correlation matrices
n9s=dblarr(n_mcmc,(nfit^2-nfit)/2+1)   ; A place to save the number of 9s in corr

; Filenames for the output
startfile=tag+'_start.dat'
truefile= tag+'_true.dat'
fitfile = tag+'_fit.dat'
errfile = tag+'_err.dat'
fomfile = tag+'_fom.dat'
corrfile= tag+'_corr.dat'
n9file= tag+'_n9.dat'

; Output space
tp=dblarr(npar_true)
fp=dblarr(npar_fit)
sp=dblarr(npar_fit)

chisq=0d
dof=0d
detcorr=0d
sigprod=0d
pegflag=0
niter=0
failflag=0

; We'll play with this a bit to figure out what the right value
; is. 2*alpha seems too small for this.  The tough thing is the
; discrepancy between the true value and the apparent size of the
; lensed image.
if not keyword_set(cutrad_scale) then cutrad_scale=2d
mincutrad=10d
maxcutrad=200d

; Set up the parameter names
fpnames=strarr(npar_fit)
tpnames=strarr(npar_true)
fpnames[0:5]=['logS0','X_c','Y_c','Alpha','e_+','e_x']
tpnames[0:5]=['logS0','X_c','Y_c','Alpha','e_+','e_x']
fpnames[npar_fit-6:npar_fit-1] = ['g_1','g_2','psi1_1','psi1_2','psi3_1','psi3_2']
tpnames[npar_true-6:npar_true-1]=['g_1','g_2','psi1_1','psi1_2','psi3_1','psi3_2']
if keyword_set(sersic_true) then tpnames[6]='n'
if keyword_set(sersic_fit) then fpnames[6]= 'n'
if keyword_set(pseudogauss_true) then tpnames[6:7]=['k1','k2']
if keyword_set(pseudogauss_fit) then fpnames[6:7]= ['k1','k2']
if keyword_set(moffat_true) then tpnames[6]='b'
if keyword_set(moffat_fit) then fpnames[6]= 'b'

; Get the true parameter ranges.
pmins=dblarr(npar_true)
pmaxs=dblarr(npar_true)

pmins[0:5]=[-1d,$                ;logS0
            -5d,-5d,$           ;X_c, Y_c
            1.5d,$              ;alpha
            -0.6d,-0.6d]        ;Eplus, Ecross

pmins[npar_true-6:npar_true-1]=[-0.5d,-0.5d,$ ;g1, g2
                                -1d-3,-1d-3,$ ;psi11,psi12
                                -1d-3,-1d-3]  ;psi31,psi32

pmaxs[0:5]=[1d,$                ;logS0
            5d,5d,$             ;X_c, Y_c
            8d,$                ;alpha
            0.6d,0.6d]          ;Eplus, Ecross

pmaxs[npar_true-6:npar_true-1]=[0.5d,0.5d,$ ;g1, g2
                                1d-3,1d-3,$ ;psi11,psi12
                                1d-3,1d-3]  ;psi31,psi32

; So that the peak is at least 1 sigma above the background level (the
; max bit is in there so that the min isn't TOO small).
if img_noise gt 0d then $
   pmins[0]=max([alog10(2d*!dpi*pmaxs[3]^2*img_noise),pmins[0]])

if keyword_set(sersic_true) then begin
; The Sersic index will run from 0.2 to 3
   pmins[6]=0.2d
   pmaxs[6]=5d
endif else if keyword_set(pseudogauss_true) then begin
; The ki parameters will run from -5 to 5
   pmins[6:7]=[0d,0d]
   pmaxs[6:7]=[5d,5d]
endif else if keyword_set(moffat_true) then begin
; The Moffat slope will run from close to 1 to 5 or so
   pmins[6]=1.2d
   pmaxs[6]=7d
endif

if keyword_set(fixcenter) then begin
   pmins[1:2]=0d
   pmaxs[1:2]=0d
endif

; Make some images and do some fits
for i=0,n_mcmc-1 do begin

   forloop_status,i,n_mcmc,flnum,label='MCMC exploring...'
 
   tp=aim_parameter_draw(pmins,pmaxs,seed,ctr_lim=pmaxs[1],e_lim=pmaxs[4],$
                        sersic=keyword_set(sersic_true),$
                        pseudogauss=keyword_set(pseudogauss_true),$
                        moffat=keyword_set(moffat_true))

   te=sqrt(total(tp[4:5]^2))
   ta=tp[3]*(1d + te)/(1d - te)
   cutrad=double(ceil(7d*ta))
   win=mk_window(ceil((2d)*cutrad))

   if not keyword_set(no_recenter_dataimg) then tp[1:2]=0d

; Make the PSF
   psf=mk_gaussian_psf(psffwhm)

; If the SIS_LENS keyword is set, then we replace the random lensmodel
; with a SIS lensmodel. This is mainly to get the COMBINATIONS of
; lensing parameters which arise in a azimuthally symmetric lens.
   if keyword_set(sis_lens) then begin

; I put the same limits on the SIS parameters by using the random
; parameters drawn above:

      psi0=sqrt(total(tp[npar_true-2:npar_true-1]^2))
      g0=sqrt(total(tp[npar_true-6:npar_true-5]^2))
      angle=atan(tp[npar_true-3],tp[npar_true-4])

      x=(1d +g0)/(2d*g0)
      theta_e=3d*g0/(4d*x*psi0)
      
      radius=x*theta_e
      
      dataimg=aim_sis_img(theta_e,[radius,angle],tp,win,psf,$
                          /polar,lensmodel=lpars,sersic=keyword_set(sersic_true))
      tp[npar_true-6:npar_true-1]=lpars
   endif else if keyword_set(nfw_lens) then begin
      pos=2.5d3*(2d*randomu(seed,2) - 1d)
      
      if (n_elements(nfw_scale) eq 0) then nfw_scale=1d

      z_lens=0.187d
      d_ang=cosmo_dist(z_lens)  ; Distance to A1689 in h_70^-1 Mpc
      
      r200=1d                   ; in Mpc - naive assumption
      c200=9d
      
      theta_s=(r200/c200)/d_ang  ; in radians
      theta_s*=(180d/!dpi)*(3600d/0.05d) ; in HST pixels

      dataimg=aim_nfw_img(theta_s,c200,z_lens,pos,tp,win,lensmodel=lpars)
      tp[npar_true-6:npar_true-1]=lpars

   endif else dataimg=aim_mk_image(tp,0d,win,psf,$
                                   sersic=keyword_set(sersic_true),$
                                   pseudogauss=keyword_set(pseudogauss_true),$
                                   moffat=keyword_set(moffat_true)) 
   dataimg+=1d-9*win ; Put a floor on the image

   bigsp=aim_start_pars(dataimg,0d,win,nsigma=0d,sigma=img_noise,seed=seed,$
                        elimit=0.9d)
   bigE=norm_vector(bigsp[4:5])
   imgscale=bigsp[3]/((1d - bigE)/(1d + bigE))
   
   
; Make the real window and dataimg
   cutrad=double(ceil(cutrad_scale*imgscale))
   cutrad=max([cutrad,mincutrad])
   cutrad=min([cutrad,maxcutrad])

   if max(size(mk_window(ceil((2d)*cutrad)),/dimensions)) lt $
      min(size(win,/dimensions)) then begin
;Only cut if the new image is smaller...
      win=mk_window(ceil((2d)*cutrad))
      ctr=bigsp[1:2] + (size(dataimg,/dimensions)-1)/2

      tp[1:2]=-bigsp[1:2]
      
      dataimg=cut_stamp(dataimg,ctr,win)
   endif

; The background level
;   bkg=0.07d*win                ; Doing things in flux units
   bkg=1d-8;*win   ; I think I can put in the win part, but idl doesn't like passing around zeros
   dataimg+=bkg

; Add in the image noise
   dataimg+=randomn(seed,size(win,/dimensions))*img_noise*win

; Do some display
   parameter_status,tp,pnames=tpnames,parnum,$
                    label='Current True Parameter Values'
   disp_scaled_image,dataimg-bkg*win,imnum,title='Data Image'
; Estimate starting parameters
   sp=aim_start_pars(dataimg,bkg,win,nsigma=1d,sigma=img_noise,seed=seed,$
                    sersic=keyword_set(sersic_fit),$
                    pseudogauss=keyword_set(pseudogauss_fit),$
                    moffat=keyword_set(moffat_fit),$
                    elimit=0.9d)


   if keyword_set(starttrueflex) then sp[npar_fit-4:npar_fit-1]=tp[npar_fit-4:npar_fit-1]

; Fit the image;;;;;;;;;;;;;;;;;;;

   fixedp=-1
   nfixed=0
   if keyword_set(fixcenter) then begin
      fixedp=[1,2]
      fixedv=[0d,0d]
      nfixed=2
   endif
; If we fix k1 or k2, set it to unity (Gaussian)
   if keyword_set(pseudogauss_fit) and keyword_set(fixk1) then begin
      fixedp=set_union(fixedp,6,count=nfixed)
      if nfixed gt 1 then fixedv=[fixedv,1d] else fixedv=1d
   endif
   if keyword_set(pseudogauss_fit) and keyword_set(fixk2) then begin
      fixedp=set_union(fixedp,7,count=nfixed)
      if nfixed gt 1 then fixedv=[fixedv,1d] else fixedv=1d
   endif
; Deal with the shearfixing
   if keyword_set(fixshear) then begin
      if abs(aim_detjacobian(tp)) gt detJaclim then begin
         sh=tp[npar_true-6:npar_true-5] + shearnoise*randomn(seed,2)
      endif else sh=[0d,0d]
      fixedp=set_union(fixedp,[npar_fit-6,npar_fit-5],count=nfixed)
      if nfixed gt 2 then fixedv=[fixedv,sh] else fixedv=sh
   endif
; Allow for a 10 parameter model with fixed psi3n=0
   if keyword_set(nopsi3) then begin
      fixedp=set_union(fixedp,[npar_fit-2,npar_fit-1],count=nfixed)
      if nfixed gt 2 then fixedv=[fixedv,0d,0d] else fixedv=[0d,0d]
   endif

   if nfixed gt 0 then sp[fixedp]=fixedv

   fp=aim_fit_image(dataimg,bkg,win,psf,$
                    startpars=sp, fixedpars=fixedp, fixedvals=fixedv,$
                    parerrors=dfp, chisq=chisq, dof=dof, corrmat=corr,$
                    detcorr=detcorr, sigsqproduct=sigprod,$
                    pegged=pegged, niter=niter,fail=failflag,$
                    sersic=keyword_set(sersic_fit),$
                    pseudogauss=keyword_set(pseudogauss_fit),$
                    moffat=keyword_set(moffat_fit),$
                    img_error=double(keyword_set(const_error))*win*img_noise,$
                    maxiter=maxiter)

   pegflag=long(total(2^pegged))
   badflag=((niter ge maxiter) or (pegflag gt 0) or (failflag gt 0))

   if ((nsave gt 0) and ((i mod nsave) eq 0)) or (badflag gt 0) then begin
      fitimg=aim_mk_image(fp,bkg,win,psf,$
                         sersic=keyword_set(sersic_fit),$
                         pseudogauss=keyword_set(pseudogauss_fit),$
                         moffat=keyword_set(moffat_fit))
      
      rimg=dataimg-fitimg
      ximg=rimg^2/fitimg
      zeroes=where(fitimg eq 0,nzeroes)
      if nzeroes gt 0 then ximg[zeroes]=0d

      maxct=floor(alog10(n_mcmc))+1
      curct=floor(alog10(max([i,1])))+1
      gap=''
      for j=curct,maxct-1 do gap+='0'

      mwrfits,dataimg,$
              tag+'_dimg'+gap+strcompress(string(long(i)),/remove_all)+'.fits'
      mwrfits,fitimg,$
              tag+'_fimg'+gap+strcompress(string(long(i)),/remove_all)+'.fits'
      mwrfits,fitimg-dataimg,$
              tag+'_rimg'+gap+strcompress(string(long(i)),/remove_all)+'.fits'
      mwrfits,ximg,$
              tag+'_ximg'+gap+strcompress(string(long(i)),/remove_all)+'.fits'

      dataimg*=0d
      fitimg*=0d
      rimg*=0d
      ximg*=0d
      win*=0d

   endif

   corr=convert_sym_matrix(corr)
   n9=double(floor(abs(alog10(1d -abs(corr))),/l64))
   
   d1=sqrt(total(fp[8:9]^2))*sp[3]
   d3=sqrt(total(fp[10:11]^2))*sp[3]
   psi1_err=sqrt(total(dfp[8:9]^2))
   psi3_err=sqrt(total(dfp[10:11]^2))

   corrs[i,*]=[i,corr]
   n9s[i,*]=[i,n9]
   startps[i,*]=[i,sp]
   trueps[i,*] =[i,tp]
   errps[i,*] = [i,dfp]
   fitps[i,*] = [i,fp]
   foms[i,*]=[i,chisq,dof,chisq/dof,$
              d1,d3,psi1_err,psi3_err,$
              niter,pegflag,failflag,badflag]


endfor


fitparlist='N, logS0, Xc, Yc, alpha, Eplus, Ecross, '
if keyword_set(sersic_fit) then fitparlist+='n, '
if keyword_set(pseudogauss_fit) then fitparlist+='k1, k2, '
if keyword_set(moffat_fit) then fitparlist+='b, '
fitparlist+='g_1, g_2, psi1_1, psi1_2, psi3_1, psi3_2'

trueparlist='N, logS0, Xc, Yc, alpha, Eplus, Ecross, '
if keyword_set(sersic_true) then trueparlist+='n, '
if keyword_set(pseudogauss_true) then trueparlist+='k1, k2, '
if keyword_set(moffat_true) then trueparlist+='b, '
trueparlist+='g_1, g_2, psi1_1, psi1_2, psi3_1, psi3_2'


spcom='Starting parameter values: '+fitparlist
tpcom='True parameter values: '+trueparlist
fpcom='Fit parameter values: '+fitparlist
pecom='Fit parameter errors: '+fitparlist

fomcom='Figures of Merit: '+$
       'N, chi^2, DoF, chi^2/DoF, D1, D3, sig(psi1), sig(psi3), '+$
       'Niter, PegFlag, FailFlag, Badflag'
corrcom='Correlation matrices: N, [corr mat]. Only unique elements are stored.'
n9com='9s in each Correlation matrix entry: '+$
      'N, [corr mat 9s]. Only unique elements are stored.'

save_data,startps,startfile,comment=spcom
save_data,trueps,truefile,comment=tpcom
save_data,fitps,fitfile,comment=fpcom
save_data,errps,errfile,comment=pecom
save_data,foms,fomfile,comment=fomcom
save_data,corrs,corrfile,comment=corrcom
save_data,n9s,n9file,comment=n9com

forloop_status,0,0,flnum,/delete
disp_scaled_image,0,imnum,/delete
parameter_status,0,parnum,/delete

aim_final_status,foms,8,9,10,maxiter=maxiter

print,'Started:  '+starttime
print,'Finished: '+systime()

end
