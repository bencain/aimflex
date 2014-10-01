pro MCMC_TEST, N_TEST, TAG=TAG, $
               SCALE=SCALE, THETA_E=THETA_E, RADIUS=RADIUS, PHI=PHI, $
               DIMS=DIMS, BKG=BKG, I0=I0, PSF_FWHM=PSF_FWMH, $
               TOL=TOL, $
               NOISE=NOISE, MCMC_LENS=MCMC_LENS, TRUESHEAR=TRUESHEAR
  
; This procedure runs through a number of cluster parameters and gives
; a 'solution cloud' for the parameters desired.

starttime=systime()

theta_e_ang=20d ; 20 arcsec theta_e
FWHM_ang=2d     ; 2 arcsec FWHM

if not keyword_set(scale) then scale=0.08d ;Arcsec per pixel

A=FWHM_ang/(scale*2.35d)

if not keyword_set(theta_e) then theta_e=theta_e_ang/scale else $
   theta_e=theta_e[0]/scale
if not keyword_set(radius) then radius=(1.5d)*theta_e else $
   radius=radius[0]

if not keyword_set(dims) then dims=[1L,1L]*floor(2d*1.7d*A) else $
   if n_elements(dims) ne 2 then dims=[dims[0],dims[0]]
if not keyword_set(bkg) then bkg=1d7 else bkg=bkg[0]
if not keyword_set(I0) then I0=0.3d else I0=I0[0]

if not keyword_set(tag) then tag='mcmc_test_' else tag=string(tag)+'_'
if not keyword_set(psf_fwhm) then psf_fwhm=0.25d else psf_fwhm=psf_fwhm[0]
if not keyword_set(tol) then tol=1d-13

; Make some global items
psf=mk_gaussian_psf(psf_fwhm)
win=mk_window(dims)


; Save the fits and the chisq
nsfitpars=dblarr(n_test,15)
wsfitpars=dblarr(n_test,15)
truepars=dblarr(n_test,12)

epsmin=0.3d
epsmax=0.7d

ximin=-!dpi/2d
ximax=!dpi/2d

ns_avecor=dblarr(12,12)
ns_avecor2=dblarr(12,12)
ws_avecor=dblarr(12,12)
ws_avecor2=dblarr(12,12)

usetrue=keyword_set(trueshear)


for i=0,n_test-1 do begin

   forloop_status,i,n_test,wnum,label='Fitting images...'

   I0_pct=0.15d
   I0_test=I0*(1d +(2d*randomu(seed)-1d)*I0_pct)

   A_pct=0.15d
   A_test=A*(1d +(2d*randomu(seed)-1d)*A_pct)

   eps=(epsmax-epsmin)*randomu(seed)+epsmin
   xi= (ximax - ximin)*randomu(seed)+ximin
 
   rc=1d-3*randomn(seed)
   ang=2d*!dpi*randomu(seed)
   xc=rc*cos(ang)
   yc=rc*sin(ang)


   if keyword_set(mcmc_lens) then begin
      phi=randomu(seed)*2d*!dpi
      lpars=sis_lensmodel([radius,phi],theta_e=theta_e,/polar)
   endif else lpars=dblarr(6)

   tp=convert_epars([I0_test,xc,yc,A_test,eps,xi,lpars],/ae_to_pol)
   truepars[i,*]=tp

   dataimg=mk_image(tp,bkg,win,psf,/floor)

   fp_ws=fit_image(dataimg,bkg,win,psf,chisq=chisq1,seed=seed,pnames=pnames,$
                   corrmat=wscor,dof=dof1,niter=nit1)

   if usetrue then sh_fix=tp[6:7] else sh_fix=[0d,0d]

   fp_ns=fit_image(dataimg,bkg,win,psf,chisq=chisq2,seed=seed,fixedpars=[6,7],$
                   fixedvals=sh_fix,corrmat=nscor,dof=dof2,niter=nit2)

   wsfitpars[i,*]=[fp_ws,chisq1,nit1,dof1]
   nsfitpars[i,*]=[fp_ns,chisq2,nit2,dof2]

   ws_avecor+=wscor
   ws_avecor2+=wscor^2
   ns_avecor+=nscor
   ns_avecor2+=nscor^2



endfor

ws_avecor/=double(n_test)
ns_avecor/=double(n_test)
ws_avecor2/=double(n_test)
ns_avecor2/=double(n_test)

forloop_status,0,0,wnum,/delete


; Save the data
save_data,truepars,tag+'truepars.dat',comment=$
          'True Parameters: I0, Xc, Yc, alpha, E+, Ex, '+$
          'g_1, g_2, G1_1, G1_2, G3_1, G3_2'
save_data,wsfitpars,tag+'wsfitpars.dat',comment=$
          'Fit Parameters (w/shear): I0, Xc, Yc, alpha, E+, Ex, '+$
          'g_1, g_2, G1_1, G1_2, G3_1, G3_2, chisq, Niter, dof'
save_data,nsfitpars,tag+'nsfitpars.dat',comment=$
          'Fit Parameters (no shear): I0, Xc, Yc, alpha, E+, Ex, '+$
          'g_1, g_2, G1_1, G1_2, G3_1, G3_2, chisq, Niter, dof'


save_data,ns_avecor,tag+'nsavecor.dat',comment=$
          'Mean Correlation (no shear): I0, Xc, Yc, alpha, E+, Ex, '+$
          'g_1, g_2, G1_1, G1_2, G3_1, G3_2'
save_data,ws_avecor,tag+'wsavecor.dat',comment=$
          'Mean Correlation (w/shear): I0, Xc, Yc, alpha, E+, Ex, '+$
          'g_1, g_2, G1_1, G1_2, G3_1, G3_2'
save_data,sqrt(ns_avecor2),tag+'nsavecorsq.dat',comment=$
          'Sigma of Correlation (no shear): I0, Xc, Yc, alpha, E+, Ex, '+$
          'g_1, g_2, G1_1, G1_2, G3_1, G3_2'
save_data,sqrt(ws_avecor2),tag+'wsavecorsq.dat',comment=$
          'Sigma of Correlation (w/shear): I0, Xc, Yc, alpha, E+, Ex, '+$
          'g_1, g_2, G1_1, G1_2, G3_1, G3_2'


mcmc_plots,tag,wsfitpars,nsfitpars,truepars,pnames

endtime=systime()

print, '*********************'
print, 'Ran '+strcompress(string(n_test),/remove_all)+' fits.'
print, 'Started: '+starttime
print, 'Ended:   '+endtime
print, '*********************'
print, '*       Done!       *'
print, '*********************'

end


