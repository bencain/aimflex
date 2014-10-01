pro mcmc_explore, n_mcmc, img_noise, tag, fixshear=fixshear, $
                  psffwhm=psffwhm,shearnoise=shearnoise, $
                  starttrueflex=starttrueflex, usenew=usenew
  
; This procedure makes a set of random input parameter values, creates
; a data image with the corresponding Gaussian noise properties and
; fits for the model parameters.

starttime=systime()

if keyword_set(shearnoise) then fixshear=1 else shearnoise=0d
if not keyword_set(psffwhm) then psffwhm=1.8d

trueps=dblarr(n_mcmc,12)
fitps=dblarr(n_mcmc,17) ; pars plus chisq, dof, chisq/dof, niter and fail flag
errps=dblarr(n_mcmc,12)
startps=dblarr(n_mcmc,12)


startfile=tag+'_mcmc_explore_start.dat'
truefile=tag+'_mcmc_explore_true.dat'
errfile=tag+'_mcmc_explore_err.dat'
fitfile=tag+'_mcmc_explore_fit.dat'
failfile=tag+'_mcmc_explore_failed.dat'
nfailfile=tag+'_mcmc_explore_notfailed.dat'

tp=dblarr(12)
fp=dblarr(12)
chisq=0d
niter=0d
dof=0d

cutrad_scale=2d
mincutrad=20d
maxcutrad=100d

pnames=['I0','X_c','Y_c','Alpha','e_+','e_x',$
        'g_1','g_2','G1_1','G1_2','G3_1','G3_2']

bkg=85d

nfail=0
failed=-1

pmins=[1d-3,$                   ;I0
       -3d,-3d,$                ;X_c, Y_c
       1d,$                     ;alpha
       -1d,-1d,$                ;Eplus, Ecross
       -1d,-1d,$                ;g1, g2
       -20d,-20d,$              ;G11,G12
       -20d,-20d]               ;G31,G32
   
pmaxs=[1d,$                     ;I0
       3d,3d,$                  ;X_c, Y_c
       50d,$                    ;alpha
       1d,1d,$                  ;Eplus, Ecross
       1d,1d,$                  ;g1, g2
       20d,20d,$                ;G11,G12
       20d,20d]                 ;G31,G32



; Make some fits
for i=0,n_mcmc-1 do begin

   forloop_status,i,n_mcmc,flnum,label='MCMC exploring...'
 
   tp=parameter_draw(pmins,pmaxs,seed,ctr_lim=3d,e_lim=0.9d)
          
;   tp[0]=randomu(seed,/double)*(1d - 1d-3)+1d-3 ; I0 from 1d-3 to 1
;   tp[1:2]=6d*randomu(seed,2,/double)-3d        ; Xc/Yc from -3 to 3
   
;   alpha_min=1d
;   alpha_max=10d
   
;   tp[3]=alpha_min+(alpha_max-alpha_min)*randomu(seed,/double)
   
;   mag=0.9d*randomu(seed,/double)
;   phase=2d*!dpi*randomu(seed,/double)
;   tp[4:5]=mag*[cos(phase),sin(phase)]

;   mag=0.9d*randomu(seed,/double)
;   phase=2d*!dpi*randomu(seed,/double)
;   tp[6:7]=mag*[cos(phase),sin(phase)]    ; g_i from -1 to 1

;   mag=20d*randomu(seed,/double)
;   phase=2d*!dpi*randomu(seed,/double)
;   tp[8:9]=mag*[cos(phase),sin(phase)]    ; G_1i from -20 to 20

;   mag=20d*randomu(seed,/double)
;   phase=2d*!dpi*randomu(seed,/double)
;   tp[10:11]=mag*[cos(phase),sin(phase)]  ; G_3i from -20 to 20

                                ; The limits on Gn correspond to
                                ; limiting each of the components to
                                ; less than 0.4 1/arcsec

   parameter_status,tp,pnames=pnames,parnum,$
                    label='Current True Parameter Values'
   
   trueps[i,*]=tp
  
   cutrad=double(ceil(cutrad_scale*tp[3]))
   cutrad=max([cutrad,mincutrad])
   cutrad=min([cutrad,maxcutrad])

   win=mk_window(ceil((2d)*cutrad))
   psf=mk_gaussian_psf(psffwhm)

   if img_noise le 0d then begin
      dataimg=mk_image(tp,bkg,win,psf)
   endif else begin
      dataimg=mk_image(tp,bkg,win,psf,/noise,gaussian_noise=img_noise)
   endelse

   disp_scaled_image,dataimg-bkg*win,imnum,title='Data Image'

   sp=start_pars(dataimg,bkg,nsigma=0.5d,sigma=img_noise*bkg,seed=seed)

   if keyword_set(starttrueflex) then sp[8:11]=tp[8:11]

; Fit the image;;;;;;;;;;;;;;;;;;;
   if keyword_Set(usenew) then begin

      if keyword_set(fixshear) then begin
         fs=tp[6:7] + shearnoise*randomn(seed,2)
         sp[6:7]=fs
         fp=new_fit_image(dataimg,bkg,win,psf,startpars=sp,chisq=chisq,dof=dof,$
                      parerrors=dfp, fixedpars=[6,7], fixedvals=fs,$
                      niter=niter,ntry=1, fail=fail)
      endif else begin
         fp=new_fit_image(dataimg,bkg,win,psf,startpars=sp,chisq=chisq,dof=dof,$
                      parerrors=dfp,niter=niter,ntry=1,fail=fail)
      endelse

   endif else begin

      if keyword_set(fixshear) then begin
         fs=tp[6:7] + shearnoise*randomn(seed,2)
         sp[6:7]=fs
         fp=fit_image(dataimg,bkg,win,psf,startpars=sp,chisq=chisq,dof=dof,$
                      parerrors=dfp, fixedpars=[6,7], fixedvals=fs,$
                      niter=niter,ntry=3, fail=fail)
      endif else begin
         fp=fit_image(dataimg,bkg,win,psf,startpars=sp,chisq=chisq,dof=dof,$
                      parerrors=dfp,niter=niter,ntry=3,fail=fail)
      endelse

endelse




   if fail then begin
      failed=set_union(failed,i)
      nfail++
   endif

   startps[i,*]=sp
   errps[i,*]=dfp
   fitps[i,*]=[fp,chisq,dof,chisq/dof,niter,fail]

endfor

failed=failed[sort(failed)]
notfailed=set_difference(lindgen(n_mcmc),failed)

scom='Starting parameter values: '+$
     'I0, Xc, Yc, alpha, Eplus, Ecross, g_1, g_2, G1_1, G1_2, G3_1, G3_2'
tcom='True parameter values: '+$
     'I0, Xc, Yc, alpha, Eplus, Ecross, g_1, g_2, G1_1, G1_2, G3_1, G3_2'
ecom='Fit parameter errors: '+$
     'I0, Xc, Yc, alpha, Eplus, Ecross, g_1, g_2, G1_1, G1_2, G3_1, G3_2'
fcom='Fit parameter values and figures of merit: '+$
     'I0, Xc, Yc, alpha, Eplus, Ecross, g_1, g_2, G1_1, G1_2, G3_1, G3_2, '+$
     'chi^2, DoF, chi^2/DoF, N_iter, Fail Flag'
failcom=$
   'Fit parameter values and figures of merit (failed fits): '+$
   'I0, Xc, Yc, alpha, Eplus, Ecross, g_1, g_2, G1_1, G1_2, G3_1, G3_2, '+$
   'chi^2, DoF, chi^2/DoF, N_iter, Fail Flag'
nfailcom=$
   'Fit parameter values and figures of merit (successful fits): '+$
   'I0, Xc, Yc, alpha, Eplus, Ecross, g_1, g_2, G1_1, G1_2, G3_1, G3_2, '+$
   'chi^2, DoF, chi^2/DoF, N_iter, Fail Flag'

save_data,startps,startfile,comment=scom
save_data,errps,errfile,comment=ecom
save_data,fitps,fitfile,comment=fcom
save_data,trueps,truefile,comment=tcom
save_data,fitps[notfailed,*],nfailfile,comment=nfailcom
if nfail gt 0 then save_data,fitps[failed,*],failfile,comment=failcom

forloop_status,0,0,flnum,/delete
disp_scaled_image,0,imnum,/delete
parameter_status,0,parnum,/delete

print,'Done!'
print,'Started'
print,starttime
print,'Finished'
print,systime()

end
