pro COV_CASE_STUDIES, A, EPS, XI, THETA_E=THETA_E, RADIUS=RADIUS, PHI=PHI,$
                      DIMS=DIMS, BKG=BKG, I0=I0, XC=XC, YC=YC, TAG=TAG, $
                      PSF_FWHM=PSF_FWMH, TOL=TOL, NOISE=NOISE, GNOISE=GNOISE,$
                      SKIPPOL=SKIPPOL, SKIPAE=SKIPAE

; This procedure creates images of size DIMS with and without lensing
; by an SIS of size THETA_E.  The image center is located at a
; distance of RADIUS and an angle PHI from the SIS center.

if not keyword_set(theta_e) then theta_e=200d else theta_e=theta_e[0]
if not keyword_set(radius) then radius=(1.5d)*theta_e else radius=radius[0]
if not keyword_set(phi) then phi=0d else phi=phi[0]
if not keyword_set(dims) then dims=[1L,1L]*floor((2d)*(1.7d)*max(A)) else $
if n_elements(dims) ne 2 then dims=[dims[0],dims[0]]
if not keyword_set(bkg) then bkg=10000000d else bkg=bkg[0]
if not keyword_set(I0) then I0=0.3d else I0=I0[0]
if not keyword_set(Xc) then Xc=0d else Xc=Xc[0]
if not keyword_set(Yc) then Yc=0d else Yc=Yc[0]
if not keyword_set(tag) then tag='case_study_' else tag=string(tag)+'_'
if not keyword_set(psf_fwhm) then psf_fwhm=0.25d else psf_fwhm=psf_fwhm[0]
if not keyword_set(tol) then tol=1d-13
if not keyword_set(gnoise) then gnoise=0d

psf=mk_gaussian_psf(psf_fwhm)
win=mk_window(dims)

; Make sure we have enough of each parameter
ngpar=min([n_elements(A),n_elements(eps),n_elements(xi)])

; Make the lens parameters
lpars=sis_lensmodel([radius,phi],theta_e=theta_e,/polar)

; Create an array for the true AE pars
ae_pars=dblarr(2*ngpar,13)
ae_pars[*,0]=dindgen(2*ngpar)+1d

ae_pars[*,1]=I0
ae_pars[*,2]=Xc
ae_pars[*,3]=Yc

ae_pars[*,4]=[A[0:ngpar-1],A[0:ngpar-1]]
ae_pars[*,5]=[eps[0:ngpar-1],eps[0:ngpar-1]]
ae_pars[*,6]=[xi[0:ngpar-1],xi[0:ngpar-1]]

ae_pars[*,7]=lpars[0]*[dblarr(ngpar),(dblarr(ngpar)+1d)]
ae_pars[*,8]=lpars[1]*[dblarr(ngpar),(dblarr(ngpar)+1d)]
ae_pars[*,9]=lpars[2]*[dblarr(ngpar),(dblarr(ngpar)+1d)]
ae_pars[*,10]=lpars[3]*[dblarr(ngpar),(dblarr(ngpar)+1d)]
ae_pars[*,11]=lpars[4]*[dblarr(ngpar),(dblarr(ngpar)+1d)]
ae_pars[*,12]=lpars[5]*[dblarr(ngpar),(dblarr(ngpar)+1d)]

; Create an array for the true POL pars
pol_pars=dblarr(2*ngpar,13)
pol_pars[*,0]=dindgen(2*ngpar)+1d
for i=0,2*ngpar-1 do pol_pars[i,1:12]=$
   convert_epars(ae_pars[i,1:12],/ae_to_pol)


; Save the true parameters and global parameters
save_data,pol_pars,tag+'pol_truepars.dat',$
          comment='POL true parameters: N, I0, Xc, Yc, alpha, e+, ex, '+$
          'g1,g2, G11, G12, G31, G32'

save_data,ae_pars,tag+'ae_truepars.dat',$
          comment='AE true parameters: N, I0, Xc, Yc, A, eps, xi, '+$
          'g1,g2, G11, G12, G31, G32'


save_data,transpose([bkg,psf_fwhm,theta_e,radius,phi]),tag+'globalpars.dat',$
          comment='Global Parameters: BKG, PSF_FWHM, THETA_E, R_SIS, PHI_SIS'

; Create storage for results and errors

pol_fits=dblarr(2*ngpar,13)
pol_fits_com='POL fit results: N, I0, Xc, Yc, alpha, e+, ex, '+$
             'g1,g2, G11, G12, G31, G32'

pol_fig_of_m=dblarr(2*ngpar,5)
pol_fig_of_m_com='POL fit figures of merit: N, chisq, dof, chisq/dof, Niter'

pol_errs=dblarr(2*ngpar,13)
pol_errs_com='POL fit errors: N, I0, Xc, Yc, alpha, e+, ex, '+$
             'g1,g2, G11, G12, G31, G32'

pol_cov_com='Covariance matrix: I0, Xc, Yc, alpha, e+, ex, '+$
            'g1,g2, G11, G12, G31, G32'
pol_corr_com='Correlation matrix: I0, Xc, Yc, alpha, e+, ex, '+$
            'g1,g2, G11, G12, G31, G32'

pol_start=dblarr(2*ngpar,13)
pol_start_com='POL starting parameters: N, I0, Xc, Yc, alpha, e+, ex, '+$
             'g1,g2, G11, G12, G31, G32'

; Create images and fit with the POL method

for i=3*ngpar*keyword_set(skippol),2*ngpar-1 do begin

   forloop_status,i,2*ngpar,wnum,label='POL fitting'

   if i+1 lt 10 then fill='00' else $
      if i+1 lt 100 then fill='0' else fill=''

   casetag=tag+'pol'+fill+strcompress(string(i+1),/remove_all)+'_'

; Make the data image and save it.
   img=mk_image(pol_pars[i,1:12],bkg,win,psf,/floor,noise=keyword_set(noise),$
                gnoise=gnoise)
   mwrfits,img,tag+'dataimg'+fill+strcompress(string(i+1),/remove_all)+'.fits'

; Fit with the POL parameter set
   niter=0
   while niter lt 10 do begin
      delvarx,startps
      fit=fit_image(img,bkg,win,psf,parerrors=err,chisq=chisq,dof=dof,$
                    covmat=cov,corrmat=corr,niter=niter,tol=tol,$
                    startpars=startps,seed=seed)
   endwhile

   fitimg=mk_image(fit,bkg,win,psf)


; Store the results
   pol_fits[i,*]=[i+1,fit]
   pol_fig_of_m[i,*]=[i+1,chisq,dof,chisq/dof,niter]
   pol_errs[i,*]=[i+1,err]
   pol_start[i,*]=[i+1,startps]

   resid=img-fitimg

   mwrfits,resid,casetag+'resid.fits'
   mwrfits,fitimg,casetag+'fitimg.fits'

   save_data,cov,casetag+'cov.dat',comment=casetag+'cov.dat '+pol_cov_com
   save_data,corr,casetag+'corr.dat',$
             comment=casetag+'corr.dat '+pol_corr_com,delimiter=' & ',$
             tail='\\', format='(f6.3)'

endfor

save_data,pol_fits,tag+'pol_fits.dat',comment=pol_fits_com
save_data,pol_fig_of_m,tag+'pol_fom.dat',comment=pol_fig_of_m_com
save_data,pol_errs,tag+'pol_errs.dat',comment=pol_errs_com
save_data,pol_start,tag+'pol_start.dat',comment=pol_start_com

; Now do the AE fitting


; Create storage for results and errors

ae_fits=dblarr(2*ngpar,13)
ae_fits_com='AE fit results: N, I0, Xc, Yc, A, eps, xi, '+$
             'g1,g2, G11, G12, G31, G32'

ae_fig_of_m=dblarr(2*ngpar,5)
ae_fig_of_m_com='AE fit figures of merit: N, chisq, dof, chisq/dof, Niter'

ae_errs=dblarr(2*ngpar,13)
ae_errs_com='AE fit errors: N, I0, Xc, Yc, A, eps, xi, '+$
             'g1,g2, G11, G12, G31, G32'

ae_cov_com='Covariance matrix: I0, Xc, Yc, A, eps, xi, '+$
            'g1,g2, G11, G12, G31, G32'
ae_corr_com='Correlation matrix: I0, Xc, Yc, A, eps, xi, '+$
            'g1,g2, G11, G12, G31, G32'

ae_start=dblarr(2*ngpar,13)
ae_start_com='AE starting parameters: N, I0, Xc, Yc, A, eps, xi, '+$
             'g1,g2, G11, G12, G31, G32'

; Create images and fit with the AE method

for i=3*ngpar*keyword_set(skipae),2*ngpar-1 do begin

   forloop_status,i,2*ngpar,wnum,label='AE fitting'

   if i+1 lt 10 then fill='00' else $
      if i+1 lt 100 then fill='0' else fill=''

   casetag=tag+'ae'+fill+strcompress(string(i+1),/remove_all)+'_'

; Make the data image (but we've already saved it).
   img=mk_image(pol_pars[i,1:12],bkg,win,psf,/floor,noise=keyword_set(noise),$
                gnoise=gnoise)

; Fit with the AE parameter set
   niter=0
   while niter lt 10 do begin
      delvarx,startps
      fit=fit_image(img,bkg,win,psf,parerrors=err,chisq=chisq,dof=dof,$
                    covmat=cov,corrmat=corr,niter=niter,startpars=startps,$
                    /ae_fit,seed=seed)
   endwhile

   fitimg=mk_image(convert_epars(fit,/ae_to_pol),bkg,win,psf)

; Store the results
   ae_fits[i,*]=[i+1,fit]
   ae_fig_of_m[i,*]=[i+1,chisq,dof,chisq/dof,niter]
   ae_errs[i,*]=[i+1,err]
   
   ae_start[i,*]=[i+1,startps]

   resid=img-fitimg

   mwrfits,resid,casetag+'resid.fits'
   mwrfits,fitimg,casetag+'fitimg.fits'

   save_data,cov,casetag+'cov.dat',comment=casetag+'cov.dat '+ae_cov_com
   save_data,corr,casetag+'corr.dat',$
             comment=casetag+'corr.dat '+ae_corr_com, delimiter=' & ',$
             tail='\\', format='(f6.3)'

endfor

forloop_status,0,0,wnum,/delete

save_data,ae_fits,tag+'ae_fits.dat',comment=ae_fits_com
save_data,ae_fig_of_m,tag+'ae_fom.dat',comment=ae_fig_of_m_com
save_data,ae_errs,tag+'ae_errs.dat',comment=ae_errs_com
save_data,ae_start,tag+'ae_start.dat',comment=ae_start_com


print,'All done!!!'

end
