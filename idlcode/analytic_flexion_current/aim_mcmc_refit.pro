pro aim_mcmc_refit, TAG, REFITS, IMG_NOISE,$
                   TAG_ADDON=TAG_ADDON,$
                   ORIG_START=ORIG_START,$
                                ; The file tag in TAG for the original
                                ; set of AIM_MCMC_EXPLORE data.
                                ;
                                ; The indices of each set of
                                ; parameters to be refit in REFITS.
                                ; If REFITS isn't specified,
                                ; then all are refit.  
                                ;
                                ; TAG_ADDON allows for an additional
                                ; retagging. 
                                ;
                                ; If ORIG_START is set then the
                                ; original starting parameters are
                                ; used as well.
                                ;
                                ; The rest of the
                                ; parameters are as in AIM_MCMC_EXPLORE
                   fixshear=fixshear, fixcenter=fixcenter,$
                   psffwhm=psffwhm,shearnoise=shearnoise, $
                   cutrad_scale=cutrad_scale, $
                   starttrueflex=starttrueflex,$
                   no_recenter_dataimg=no_recenter_dataimg
;;;;;;;;

orig_startps=read_data(tag+'_start.dat',size=dims,/quiet)
orig_trueps=read_data(tag+'_true.dat',/quiet)

if ((n_elements(refits) eq 0) or $
    (max(refits) ge dims[0]) or $
    (min(refits) lt 0)) then refits=lindgen(dims[0])

n_mcmc=n_elements(refits)

if not keyword_set(tag_addon) then tag_addon=''
retag=tag+'_refit'+tag_addon

starttime=systime()

if keyword_set(shearnoise) then fixshear=1 else shearnoise=0d
if not keyword_set(psffwhm) then psffwhm=1.8d

startps=dblarr(n_mcmc,13)       ; An index number plus starting parameters
trueps=dblarr(n_mcmc,13)        ; An index number plus true parameters
fitps=dblarr(n_mcmc,13)         ; An index number plus fit parameters
errps=dblarr(n_mcmc,13)         ; An index number plus parameter errors
foms=dblarr(n_mcmc,9)           ; An index number plus figures of merit:
                                ;   chisq, DoF, chisq/DoF, detCORR,
                                ;   prod(sigsq), npegged, niter, failflag

; Filenames for the output
startfile=retag+'_start.dat'
truefile= retag+'_true.dat'
fitfile = retag+'_fit.dat'
convtruefile= retag+'_convtrue.dat'
convfitfile = retag+'_convfit.dat'
errfile = retag+'_err.dat'
fomfile = retag+'_fom.dat'

; Output space
tp=dblarr(12)
fp=dblarr(12)
sp=dblarr(12)

chisq=0d
dof=0d
detcorr=0d
sigprod=0d
npeg=0
niter=0
failflag=0

; We'll play with this a bit to figure out what the right value
; is. 2*alpha seems too small for this.  The tough thing is the
; discrepancy between the true value and the apparent size of the
; lensed image.
if not keyword_set(cutrad_scale) then cutrad_scale=2d
mincutrad=20d
maxcutrad=200d

pnames=['logN0','X_c','Y_c','Alpha','e_+','e_x',$
        'g_1','g_2','psi1_1','psi1_2','psi3_1','psi3_2']
convpnames=['logN0','X_c','Y_c','Alpha','e_+','e_x',$
            'g_1','g_2','G1_1','G1_2','G3_1','G3_2']


bkg=85d

pmins=[2d,$                   ;logN0
       -5d,-5d,$              ;X_c, Y_c
       1d,$                   ;alpha
       -1d,-1d,$              ;Eplus, Ecross
       -0.9d,-0.9d,$          ;g1, g2
       -2d-2,-2d-2,$          ;psi11,psi12
       -2d-2,-2d-2]           ;psi31,psi32
   
pmaxs=[6d,$                   ;logN0
       5d,5d,$                ;X_c, Y_c
       50d,$                  ;alpha
       0.8d,0.8d,$            ;Eplus, Ecross
       0.9d,0.9d,$            ;g1, g2
       2d-2,2d-2,$            ;psi11,psi12
       2d-2,2d-2]             ;psi31,psi32

if keyword_set(fixcenter) then begin
   pmins[1:2]=0d
   pmaxs[1:2]=0d
endif

; Make some images and do some fits
for i=0,n_mcmc-1 do begin

   forloop_status,i,n_mcmc,flnum,label='Refitting...'

   ; Put in the true parameters from the original fit
   tp=transpose(orig_trueps[refits[i],1:12])

   parameter_status,tp,pnames=pnames,parnum,$
                    label='Current True Parameter Values'

   cutrad=double(ceil(cutrad_scale*tp[3]))
   cutrad=max([cutrad,mincutrad])
   cutrad=min([cutrad,maxcutrad])

   win=mk_window(ceil((2d)*cutrad))
   psf=mk_gaussian_psf(psffwhm)


   if img_noise le 0d then begin
      dataimg=aim_mk_image(tp,bkg,win,psf)
   endif else begin
      dataimg=aim_mk_image(tp,bkg,win,psf,gaussian_noise=img_noise)
   endelse

      

   disp_scaled_image,dataimg-bkg*win,imnum,title='Data Image'
   
; Put in the starting parameters
   if keyword_set(orig_start) then begin
      sp=transpose(orig_startps[refits[i],1:12])
   endif else begin
      sp=aim_start_pars(dataimg,bkg,win,nsigma=0.5d,sigma=img_noise*bkg,seed=seed)
   endelse

   if keyword_set(starttrueflex) then sp[8:11]=tp[8:11]

; Fit the image;;;;;;;;;;;;;;;;;;;

   if keyword_set(fixshear) then begin
      fixedv=tp[6:7] + shearnoise*randomn(seed,2)
      sp[6:7]=fixedv
      fixedp=[6,7]
      if keyword_set(fixcenter) then begin
         fixedp=[1,2,fixedp]
         fixedv=[0d,0d,fixedv]
      endif

      fp=aim_fit_image(dataimg,bkg,win,psf,$
                      startpars=sp, fixedpars=fixedp, fixedvals=fixedv,$
                      parerrors=dfp, chisq=chisq, dof=dof,$
                      detcorr=detcorr, sigsqproduct=sigprod,$
                      npegged=npeg, niter=niter,fail=failflag)
   endif else begin
      fp=aim_fit_image(dataimg,bkg,win,psf,$
                       startpars=sp, $
                       parerrors=dfp, chisq=chisq, dof=dof,$
                       detcorr=detcorr, sigsqproduct=sigprod,$
                       npegged=npeg, niter=niter,fail=failflag)
   endelse

   startps[i,*]=[refits[i],sp]
   trueps[i,*] =[refits[i],tp]
   errps[i,*] = [refits[i],dfp]
   fitps[i,*] = [refits[i],fp]
   foms[i,*]=[refits[i],chisq,dof,chisq/dof,detcorr,sigprod,npeg,niter,failflag]

endfor

convtps=trueps
convfps=fitps
for i=0,n_mcmc-1 do begin
   convtps[i,1:12]=convert_flexpars(convtps[i,1:12],/psin_to_gn)
   convfps[i,1:12]=convert_flexpars(convfps[i,1:12],/psin_to_gn)
endfor

spcom='Starting parameter values: N, logN0, Xc, Yc, alpha, Eplus, Ecross, '+$
      'g_1, g_2, psi1_1, psi1_2, psi3_1, psi3_2'
tpcom='True parameter values: N, logN0, Xc, Yc, alpha, Eplus, Ecross, '+$
      'g_1, g_2, psi1_1, psi1_2, psi3_1, psi3_2'
fpcom='Fit parameter values: N, logN0, Xc, Yc, alpha, Eplus, Ecross, '+$
      'g_1, g_2, psi1_1, psi1_2, psi3_1, psi3_2'
convtpcom='Converted True parameter values: '+$
          'N, logN0, Xc, Yc, alpha, Eplus, Ecross, '+$
          'g_1, g_2, G1_1, G1_2, G3_1, G3_2'
convfpcom='Converted Fit parameter values: '+$
          'N, logN0, Xc, Yc, alpha, Eplus, Ecross, '+$
          'g_1, g_2, G1_1, G1_2, G3_1, G3_2'
pecom='Fit parameter errors: N, logN0, Xc, Yc, alpha, Eplus, Ecross, '+$
      'g_1, g_2, psi1_1, psi1_2, psi3_1, psi3_2'
fomcom='Figures of Merit: '+$
       'N, chi^2, DoF, chi^2/DoF, det(CORR), prod(err^2), Npegged, Niter, FailFlag'

save_data,startps,startfile,comment=spcom
save_data,trueps,truefile,comment=tpcom
save_data,fitps,fitfile,comment=fpcom
save_data,convtps,convtruefile,comment=convtpcom
save_data,convfps,convfitfile,comment=convfpcom
save_data,errps,errfile,comment=pecom
save_data,foms,fomfile,comment=fomcom

forloop_status,0,0,flnum,/delete
disp_scaled_image,0,imnum,/delete
parameter_status,0,parnum,/delete

print,'Done!'
print,'Started:  '+starttime
print,'Finished: '+systime()

end
