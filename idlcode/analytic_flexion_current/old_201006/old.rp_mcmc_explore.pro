pro rp_mcmc_explore, n_mcmc, img_noise, tag, $
                     fixshear=fixshear, fixcenter=fixcenter,$
                     psffwhm=psffwhm,shearnoise=shearnoise, $
                     cutrad_scale=cutrad_scale, $
                     starttrueflex=starttrueflex,$
                     sis_lens=sis_lens,$
                     no_recenter_dataimg=no_recenter_dataimg,$
                     sersic=sersic

; This procedure makes a set of random input parameter values, creates
; a data image with the corresponding Gaussian noise properties and
; fits for the model parameters.

starttime=systime()

if keyword_set(shearnoise) then fixshear=1 else shearnoise=0d
if not keyword_set(psffwhm) then psffwhm=1.8d
if keyword_set(sersic) then npar=13 else npar=12

startps=dblarr(n_mcmc,npar+1)   ; An index number plus starting parameters
trueps=dblarr(n_mcmc,npar+1)    ; An index number plus true parameters
fitps=dblarr(n_mcmc,npar+1)     ; An index number plus fit parameters
errps=dblarr(n_mcmc,npar+1)     ; An index number plus parameter errors
foms=dblarr(n_mcmc,10)          ; An index number plus figures of merit:
                                ;   chisq, DoF, chisq/DoF, detCORR,
                                ;   N_corr 9s, prod(sigsq), pegflag, niter,
                                ;   failflag
nfit=npar
if keyword_set(fixshear) then nfit-=2
if keyword_set(fixcenter) then nfit-=2
corrs=dblarr(n_mcmc,(nfit^2-nfit)/2+1) ; A place to save the correlation matrices
n9s=dblarr(n_mcmc,(nfit^2-nfit)/2+1)   ; A place to save the number of 9s in corr

; Filenames for the output
startfile=tag+'_start.dat'
truefile= tag+'_true.dat'
fitfile = tag+'_fit.dat'
convtruefile= tag+'_convtrue.dat'
convfitfile = tag+'_convfit.dat'
errfile = tag+'_err.dat'
fomfile = tag+'_fom.dat'
corrfile= tag+'_corr.dat'
n9file= tag+'_n9.dat'

; Output space
tp=dblarr(npar)
fp=dblarr(npar)
sp=dblarr(npar)

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
if not keyword_set(cutrad_scale) then cutrad_scale=3d
mincutrad=20d
maxcutrad=200d

pnames=['logN0','X_c','Y_c','Alpha','e_+','e_x',$
        'g_1','g_2','psi1_1','psi1_2','psi3_1','psi3_2']
convpnames=['logN0','X_c','Y_c','Alpha','e_+','e_x',$
            'g_1','g_2','G1_1','G1_2','G3_1','G3_2']
if keyword_set(sersic) then begin
   pnames=[pnames[0:5],'n',pnames[6:11]]
   pnames=[convpnames[0:5],'n',convpnames[6:11]]
endif

bkg=85d

pmins=[2d,$                   ;logN0
       -5d,-5d,$              ;X_c, Y_c
       2d,$                  ;alpha
       -0.8d,-0.8d,$          ;Eplus, Ecross
       -0.8d,-0.8d,$          ;g1, g2
       -2d-2,-2d-2,$          ;psi11,psi12
       -2d-2,-2d-2]           ;psi31,psi32
   
pmaxs=[6d,$                   ;logN0
       5d,5d,$                ;X_c, Y_c
       50d,$                  ;alpha
       0.8d,0.8d,$            ;Eplus, Ecross
       0.8d,0.8d,$            ;g1, g2
       2d-2,2d-2,$            ;psi11,psi12
       2d-2,2d-2]             ;psi31,psi32

if keyword_set(sersic) then begin
; The Sersic index will run from 0.25 to 4
   pmins=[pmins[0:5],0.25d,pmins[6:11]]
   pmaxs=[pmaxs[0:5],4d,pmaxs[6:11]]
endif

if keyword_set(fixcenter) then begin
   pmins[1:2]=0d
   pmaxs[1:2]=0d
endif

; Make some images and do some fits
for i=0,n_mcmc-1 do begin

   forloop_status,i,n_mcmc,flnum,label='MCMC exploring...'
 
   tp=rp_parameter_draw(pmins,pmaxs,seed,ctr_lim=pmaxs[1],e_lim=pmaxs[4],$
                        sersic=keyword_set(sersic))
   if not keyword_set(no_recenter_dataimg) then tp[1:2]=0d

; If the SIS_LENS keyword is set, then we replace the random lensmodel
; with a SIS lensmodel. This is mainly to get the COMBINATIONS of
; lensing parameters which arise in a azimuthally symmetric lens.
   if keyword_set(sis_lens) then begin
                                ; Einstein radius of 10 - 100 arcsec at
                                ;                        HST resolution
      theta_e=(90d*randomu(seed)+10d)/0.05d
                                ; Image located between 1.5 and 2.5 theta_e
      radius=(randomu(seed)+1.2d)*theta_e
      angle=2d*!dpi*randomu(seed)
      lpars=sis_lensmodel([radius,angle],/polar,/psipars,$
                          theta_e=theta_e,flexionscale=1d)
      tp[npar-6:npar-1]=lpars
   endif

   cutrad=double(ceil(cutrad_scale*tp[3]))
   cutrad=max([cutrad,mincutrad])
   cutrad=min([cutrad,maxcutrad])

   win=mk_window(ceil((2d)*cutrad))
   psf=mk_gaussian_psf(psffwhm)
   
; Default will be to recenter any data image on the apparent centroid,
; but you can escape this with the keyword NO_RECENTER_DATAIMG
   if not keyword_set(no_recenter_dataimg) then begin
      if img_noise le 0d then begin
         dataimg=rp_mk_dataimage(tp,bkg,win,psf,/replace_center,sersic=keyword_set(sersic))
      endif else begin
         dataimg=rp_mk_dataimage(tp,bkg,win,psf,gaussian_noise=img_noise,$
                                 /replace_center,sersic=keyword_set(sersic))
      endelse
   endif else begin
      if img_noise le 0d then begin
         dataimg=rp_mk_image(tp,bkg,win,psf,sersic=keyword_set(sersic))
      endif else begin
         dataimg=rp_mk_image(tp,bkg,win,psf,gaussian_noise=img_noise,$
                             sersic=keyword_set(sersic))
      endelse
   endelse


; Do some display
   parameter_status,tp,pnames=pnames,parnum,$
                    label='Current True Parameter Values'
   disp_scaled_image,dataimg-bkg*win,imnum,title='Data Image'
   sp=rp_start_pars(dataimg,bkg,nsigma=0.5d,sigma=img_noise*bkg,seed=seed,$
                    sersic=keyword_set(sersic),elimit=0.9d)

   if keyword_set(starttrueflex) then sp[npar-4:npar-1]=tp[npar-4:npar-1]

; Fit the image;;;;;;;;;;;;;;;;;;;

   if keyword_set(fixshear) then begin
      fixedv=tp[npar-6:npar-5] + shearnoise*randomn(seed,2)
      sp[npar-6:npar-5]=fixedv
      fixedp=[npar-6,npar-5]
      if keyword_set(fixcenter) then begin
         fixedp=[1,2,fixedp]
         fixedv=[0d,0d,fixedv]
      endif

      fp=rp_fit_image(dataimg,bkg,win,psf,$
                      startpars=sp, fixedpars=fixedp, fixedvals=fixedv,$
                      parerrors=dfp, chisq=chisq, dof=dof, corrmat=corr,$
                      detcorr=detcorr, sigsqproduct=sigprod,$
                      pegged=pegged, niter=niter,fail=failflag,$
                      sersic=keyword_set(sersic))
   endif else begin
      fp=rp_fit_image(dataimg,bkg,win,psf,$
                      startpars=sp, $
                      parerrors=dfp, chisq=chisq, dof=dof, corrmat=corr,$
                      detcorr=detcorr, sigsqproduct=sigprod,$
                      pegged=pegged, niter=niter,fail=failflag,$
                      sersic=keyword_set(sersic))
   endelse
   
   pegflag=long(total(2^pegged))

   corr=convert_sym_matrix(corr)
   n9=double(floor(abs(alog10(1d -abs(corr))),/l64))
   corrs[i,*]=[i,corr]
   n9s[i,*]=[i,n9]
   startps[i,*]=[i,sp]
   trueps[i,*] =[i,tp]
   errps[i,*] = [i,dfp]
   fitps[i,*] = [i,fp]
   foms[i,*]=[i,chisq,dof,chisq/dof,detcorr,total(n9,/double),$
              sigprod,niter,pegflag,failflag]

endfor

convtps=trueps
convfps=fitps
for i=0,n_mcmc-1 do begin
   convtps[i,1:npar]=convert_flexpars(convtps[i,1:npar],/psin_to_gn)
   convfps[i,1:npar]=convert_flexpars(convfps[i,1:npar],/psin_to_gn)
endfor

spcom='Starting parameter values: N, logN0, Xc, Yc, alpha, Eplus, Ecross, '
if keyword_set(sersic) then spcom+='n, '
spcom+='g_1, g_2, psi1_1, psi1_2, psi3_1, psi3_2'

tpcom='True parameter values: N, logN0, Xc, Yc, alpha, Eplus, Ecross, '
if keyword_set(sersic) then tpcom+='n, '
tpcom+='g_1, g_2, psi1_1, psi1_2, psi3_1, psi3_2'

fpcom='Fit parameter values: N, logN0, Xc, Yc, alpha, Eplus, Ecross, '
if keyword_set(sersic) then fpcom+='n, '
fpcom+='g_1, g_2, psi1_1, psi1_2, psi3_1, psi3_2'

convtpcom='Converted True parameter values: '+$
          'N, logN0, Xc, Yc, alpha, Eplus, Ecross, '
if keyword_set(sersic) then convtpcom+='n, '
convtpcom+='g_1, g_2, G1_1, G1_2, G3_1, G3_2'

convfpcom='Converted Fit parameter values: '+$
          'N, logN0, Xc, Yc, alpha, Eplus, Ecross, '
if keyword_set(sersic) then convfpcom+='n, '
convfpcom+='g_1, g_2, G1_1, G1_2, G3_1, G3_2'

pecom='Fit parameter errors: N, logN0, Xc, Yc, alpha, Eplus, Ecross, '
if keyword_set(sersic) then pecom+='n, '
pecom+='g_1, g_2, psi1_1, psi1_2, psi3_1, psi3_2'

fomcom='Figures of Merit: '+$
       'N, chi^2, DoF, chi^2/DoF, det(CORR), N_corr 9s, '+$
       'prod(err^2), Niter, PegFlag, FailFlag'
corrcom='Correlation matrices: N, [corr mat]. Only unique elements are stored.'
n9com='9s in each Correlation matrix entry: '+$
      'N, [corr mat 9s]. Only unique elements are stored.'

save_data,startps,startfile,comment=spcom
save_data,trueps,truefile,comment=tpcom
save_data,fitps,fitfile,comment=fpcom
save_data,convtps,convtruefile,comment=convtpcom
save_data,convfps,convfitfile,comment=convfpcom
save_data,errps,errfile,comment=pecom
save_data,foms,fomfile,comment=fomcom
save_data,corrs,corrfile,comment=corrcom
save_data,n9s,n9file,comment=n9com

forloop_status,0,0,flnum,/delete
disp_scaled_image,0,imnum,/delete
parameter_status,0,parnum,/delete

print,'Done!'
print,'Failed fits: '+strcompress(string(n_elements(where(foms[*,9] gt 0))),/remove_all)
print,'Pegged fits: '+strcompress(string(n_elements(where(foms[*,8] gt 0))),/remove_all)
print,'Started:  '+starttime
print,'Finished: '+systime()

end
