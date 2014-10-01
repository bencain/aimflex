pro aim_fit_catalog, imgfile, red_catfile, blue_catfile, $
                     img_error=img_error, maxiter=maxiter,$
                     psf_fwhm=psf_fwhm, hst_psf=hst_psf,$
                     tag=tag, nsave=nsave,$
                     save_ximg=save_ximg, save_dimg=save_dimg,$
                     redname=redname, bluename=bluename,$
                     testlimit=testlimit, nfw_scale=nfw_scale, nfw_ctr=nfw_ctr,$
                     pseudogauss=pseudogauss, fixk1=fixk1, fixk2=fixk2, $
                     sersic=sersic, moffat=moffat,$
                     startflex=startflex, detJacmin=detJacmin, detJacmax=detJacmax,$
                     freeshear=freeshear, fit_all=fit_all, $
                     extra_shearnoise=extra_shearnoise,$ ; For use with AIM_MK_CLUSTER_IMG data
                     fit_nolens=fit_nolens,$
                     cutrad_scale=cutrad_scale,min_cutrad=min_cutrad,max_cutrad=max_cutrad

; This function takes in the catalog file output of the
; SEX_CAT_CONVERT script, the full fits image, and an output filename
; for the fit catalog.  The fitting procedure is as follows:
;
;  1. Color/Magnitude/Size discrimination
;       - Fit galaxies, member galaxies, stars, bad objects
;  2. PSF estimation
;  3. Image Fitting

; Read in the catalogs.  The entries are:
;       N, X_field, Y_field, logN0, Alpha, eplus, ecross, Bkg, Mag,
;       MagErr, Class_star  
  redcat=read_data(red_catfile,comment=rcatcom,/quiet)
  bluecat=read_data(blue_catfile,comment=bcatcom,/quiet)

  if not keyword_set(img_error) then img_error=0d
  if not keyword_set(extra_shearnoise) then extra_shearnoise=0d

; Read in the data image and make sure that it's sane
  full_image=mrdfits(imgfile)
  negpix=where(full_image le 0d,n_negpix,complement=pospix)
  if n_negpix gt 0d then full_image[negpix]=0d

; Set up lists of the fit parameters and the fixed parameters
  if (keyword_set(sersic) or keyword_set(moffat)) then n_par=13 else $
     if keyword_set(pseudogauss) then n_par=14 else n_par=12
  fitps=lindgen(n_par)
  fixps=-1
  n_fitps=n_par
  n_fixps=0
  if not keyword_set(freeshear) then begin
     fitps=set_difference(fitps,[n_par-6,n_par-5],count=n_fitps)
     fixps=set_union(fixps,[n_par-6,n_par-5],count=n_fixps)
  endif
  if keyword_set(fixk1) then begin
     fitps=set_difference(fitps,6,count=n_fitps)
     fixps=set_union(fixps,6,count=n_fixps)
  endif
  if keyword_set(fixk2) then begin
     fitps=set_difference(fitps,7,count=n_fitps)
     fixps=set_union(fixps,7,count=n_fixps)
  endif
  if keyword_set(fit_nolens) then begin
     fitps=set_difference(fitps,[1,2,lindgen(6)+n_par-6],count=n_fitps)
     fixps=set_union(fixps,[1,2,lindgen(6)+n_par-6],count=n_fixps)
     nfw_scale=0d
  endif

  if n_fixps gt 0 then fixedvals=dblarr(n_fixps)

; Pull out the magnitudes and backgrounds.  I should put zeropoints in
; at some point here.
  rmag=redcat[*,8]
  ralpha=redcat[*,4]
  rEmag=sqrt(total(redcat[*,5:6]^2,2))
  reps=sqrt((1d - rEmag)/(1d + rEmag))
  rsize=ralpha/sqrt(reps)
  bmag=bluecat[*,8]
  color=bmag-rmag
  bkg=redcat[*,7]
  n_total=n_elements(color)


; ---------- SIZE/COLOR SELECTION ----------

; We'll select galaxies to fit by size and magnitude.
;
;  fitgal=interactive_select(ralpha,rmag,aname='ALPHA',bname=redname,$
;                            prompt='Select on size/magnitude',$
;                            n_sel=n_fit,not_sel=notfit,n_not_sel=n_notfit)

  if keyword_set(fit_all) then begin
     fitgal=lindgen(n_total) 
     n_fit=n_total
  endif else fitgal=aim_cmr_select(color,rmag,rsize,n_fit_obj=n_fit)


; ---------- PSF MODEL ESTIMATION ----------

; Do something more sophisticated later.  For now, the HST PSF is
; approximately 0.09", or 1.8 pixels in FWHM.  We'll use this
; for all the galaxies, though later we'll do something more
; individualized. 

  if not keyword_set(psf_fwhm) then psf_fwhm=1.8d
  psf_model=mk_gaussian_psf(psf_fwhm)

  if keyword_set(hst_psf) then $
     psf_model=mrdfits('~/idl/data/a1689_input/psf/hst_ave_psf_rot25_small.fits')

; ---------- IMAGE FITTING ----------

; Make the data saving arrays.
  startpars=dblarr(n_fit,n_par+1) ; Index and each starting parameter.
  parfits = dblarr(n_fit,n_par+3) ; Index, object location, and each of the 
                                ; parameter fit values.
  parerrors=dblarr(n_fit,n_par+1) ; Index and each parameter error.
  meritfigs=dblarr(n_fit,14)       ; Index and the following figures of merit:
                                ;    chi^2=sum( (D-M)^2/M )
                                ;    DoF (degrees of freedom)
                                ;    chi^2/DoF
                                ;    D1=alpha*sqrt(psi11^2 + psi12^2)
                                ;    D3=alpha*sqrt(psi31^2 + psi32^2)
                                ;    sig_psi1 (error in psi11-psi12
                                ;            space)
                                ;    sig_psi3 (error in psi31-psi32
                                ;            space)
                                ;    rmag - red filter magnitude
                                ;    color = bmag - rmag
                                ;    N_iter (number of iterations to
                                ;            the returned fit)
                                ;    PegFlag (A binary representation
                                ;            of the pegged
                                ;            parameters, if any)
                                ;    FailFlag (binary fail marker)
                                ;    BadFlag (binary marker for
                                ;            failure in iterationns,
                                ;            pegging or failing.
  n2d=n_fitps
  corrmats=dblarr(n_fit,(n2d^2-n2d)/2 + 1) 
                                ; The number of parameters -2 assumes
                                ; that the shears are fixed, and if k1
                                ; and/or k2 are fixed, subtract one
                                ; for each, since CORRMAT is returned
                                ; from AIM_FIT_IMAGE as an n2d x n2d
                                ; matrix, where n2d is the number of
                                ; free parameters.

; Set up file names and their save comments
  startpars_file=tag+'_start.dat'
  parfits_file = tag+'_fit.dat'
  parerrors_file=tag+'_err.dat'
  meritfigs_file=tag+'_fom.dat'
  corrmats_file =tag+'_corr.dat'
  cfsp_parfits_file=tag+'_cfspfit.dat'

  fitparlist='logN0, Xc, Yc, alpha, Eplus, Ecross, '
  if keyword_set(sersic) then fitparlist+='n, '
  if keyword_set(pseudogauss) then fitparlist+='k1, k2, '
  if keyword_set(moffat) then fitparlist+='b, '
  
  fsp_fitparlist=fitparlist+'g_0, psi1_1, psi1_2, psi3_1, psi3_2'
  fitparlist+='g_1, g_2, psi1_1, psi1_2, psi3_1, psi3_2'
  
  startpars_com='Starting parameter values: N, '+fitparlist

  parfits_com = 'Fit parameter values: N, X_field, Y_field, '+fitparlist

  parerrors_com='Fit parameter errors: N, '+fitparlist

  meritfigs_com='Figures of Merit: '+$
                'N, chisq, DoF, chisq/DoF, D1, D3, sig_psi1, sig_psi3, '+$
                'Mag, Color, N_iter, PegFlag, FailFlag, BadFlag'

  corrmats_com='Correlation Matrices: '+$
               'N, [Converted 1D correlation matrices]'

; Some general fitting parameters
  if not keyword_set(cutrad_scale) then cutrad_scale=2d
                                ; Cut with a window with radius 
                                ; 2*alpha_init (from SExtractor)

  if not keyword_set(min_cutrad) then min_cutrad=15d  ; Minimum radius for fitting
  if not keyword_set(max_cutrad) then max_cutrad=200d ; Maximum radius for fitting

                                ; Get rid of any places where galaxies
                                ; are off the edge
  fitgal=set_intersection(fitgal,where(bkg gt 0d),count=n_fit)

  maxindex=max(redcat[fitgal,0]) ; Maximum index in the fit catalog

  if not keyword_set(nsave) then nsave=-1
  if not keyword_set(maxiter) then maxiter=1000

; Set up the shear lens model
  if not keyword_set(nfw_scale) then nfw_scale=1d
  if not keyword_set(nfw_ctr) then nfw_ctr=[3130d,2865d] ; Center of the lens in the image

  z_lens=0.183d
  d_ang=cosmo_dist(0.183)       ; Distance to A1689 in h_70^-1 Mpc

  M200=1.31d15*(0.7)^(-1)       ; From Peng et al. 2010
  c200=9.9d
  H0=70d                        ; km/s/Mpc
  Eofz=sqrt(0.3d*(1d + z_lens)^3 + 0.7d)
  Gconst=4.301d-9               ; (km/s)^2*Mpc/M_sun
  r200=(0.01d*Gconst*M200/(H0*Eofz)^2)^(1d/3d) ; Mpc

  r_s=r200/c200

  theta_s=(r_s/d_ang)*(180d/!dpi)*(3600d/0.05d) ; in HST pix

  a1689_lensmodel=nfw_lensmodel(redcat[*,1:2],z_lens,c200,theta_s,$
                                ctr=nfw_ctr,force_scale=nfw_scale)

  if not keyword_set(detJacmin) then detJacmin=0.2d
  if not keyword_set(detJacmax) then detJacmax=10d
                                ; Limit on when to not use the shear
                                ; model and instead fix the shear to
                                ; zero.  This decision is based on the
                                ; determinant of the lensing Jacobian.


; -------- Cut and fit each image ------------

  if keyword_set(testlimit) then n_fit=10
  starttime=systime()
  print,'Fitting Started: '+starttime

  clear_windows

  for i=0,n_fit-1 do begin
     if (i mod 10) eq 0 then $
        forloop_status,i,n_fit,countwnum,label='-=< Fitting the selected galaxies >=-'
   
; Create the window
     se_alpha=redcat[fitgal[i],4]
     se_Emag= sqrt(total(redcat[fitgal[i],5:6]^2))
     se_eps = sqrt((1d - se_Emag)/(1d + se_Emag))
     
;     cutrad=cutrad_scale*redcat[fitgal[i],4]
     cutrad=cutrad_scale*se_alpha/sqrt(se_eps)
     cutrad=max([cutrad,min_cutrad])
     cutrad=min([cutrad,max_cutrad])

     win=mk_window(2d*cutrad)

; Build the PSF from a model
; For now, that model is just a uniform Gaussian
     psf=psf_model              ;Eventually this will be a function

; Cut and display the stamp
     dimg=cut_stamp(full_image,redcat[fitgal[i],1:2],win)
     disp_scaled_image,alog(dimg-bkg[fitgal[i]]*win + 1d),$
                       imgwnum,title='Current Image',label='Object #'+$
                       strcompress(string(long(redcat[fitgal[i],0])),/remove_all)

; Decide about whether to use the shear model or to just set shear to zero.
     if (not keyword_set(freeshear)) then begin
        if not keyword_set(fit_nolens) then begin
           dJ=abs(aim_detjacobian(a1689_lensmodel[fitgal[i],*],scale=ralpha[fitgal[i]]))
           if (dJ gt detJacmin) and (dJ lt detJacmax) then begin
              fixedvals[where(fixps eq n_par-6)]=a1689_lensmodel[fitgal[i],0] + $
                                         extra_shearnoise*randomn(seed)
              fixedvals[where(fixps eq n_par-5)]=a1689_lensmodel[fitgal[i],1] + $
                                         extra_shearnoise*randomn(seed)
           endif else begin
              fixedvals[where(fixps eq n_par-6)]=0d
              fixedvals[where(fixps eq n_par-5)]=0d
           endelse
        endif  
     endif


; Set up parameter fixing for the Pseudo-Gaussian model
     if (keyword_set(fixk2) and keyword_set(pseudogauss)) then begin
        fixedvals[where(fixps eq 7)]=1d
     endif
     if (keyword_set(fixk1) and keyword_set(pseudogauss)) then begin
        fixedvals[where(fixps eq 6)]=1d
     endif

; If we're inputting a shear, then we can input a starting
; flexion consistent with the lens model.
;   if (keyword_set(startflex) and $
;       (sqrt(total(a1689_lensmodel[fitgal[i],0:1]^2)) lt shearlim)) then $
;          sp[n_par-4:n_par-1]=a1689_lensmodel[fitgal[i],2:5]

     delvarx,sp

     fp=aim_fit_image(dimg,bkg[fitgal[i]],win,psf,$
                     img_error=img_error, maxiter=maxiter,$
                     startpars=sp,fixedpars=fixps,fixedvals=fixedvals,$
                     parerrors=pe,chisq=chisq,dof=dof, seed=seed,$
                     pegflag=pegflag, niter=niter,failflag=failflag,$
                     corrmat=corr,$
                     sersic=keyword_set(sersic), $
                     pseudogauss=keyword_set(pseudogauss),$
                     moffat=keyword_set(moffat),$
                     fit_index=redcat[fitgal[i],0])

     badflag=((pegflag gt 0) or (niter ge maxiter) or (failflag gt 0))

; Save the data, fit and residual images if i is a multiple of nsave.
     if ((nsave gt 0) and ((i mod nsave) eq 0)) or (badflag ne 0) then begin
        fimg=aim_mk_image(fp,bkg[fitgal[i]],win,psf,$
                          sersic=keyword_set(sersic),$
                          pseudogauss=keyword_set(pseudogauss),$
                          moffat=keyword_set(moffat))
        nans=where(finite(fimg) eq 0,nnans)
        zeroes=where(fimg eq 0d,nzeroes)
        if nnans gt 0d then fimg[nans]=9d9
        rimg=fimg-dimg
        ximg=rimg^2/fimg
        if nzeroes gt 0d then ximg[zeroes]=0d
        maxct=floor(alog10(maxindex))+1
        curct=floor(alog10(redcat[fitgal[i],0]))+1
        gap=''
        for j=curct,maxct-1 do gap+='0'
        
        mwrfits,fimg,tag+'_fimg'+gap+$
                strcompress(string(long(redcat[fitgal[i],0])),/remove_all)+'.fits'
        mwrfits,rimg,tag+'_rimg'+gap+$
                strcompress(string(long(redcat[fitgal[i],0])),/remove_all)+'.fits'

        if keyword_set(save_dimg) then $
           mwrfits,dimg,tag+'_dimg'+gap+$
                   strcompress(string(long(redcat[fitgal[i],0])),/remove_all)+'.fits'
        if keyword_set(save_ximg) then $
           mwrfits,ximg,tag+'_ximg'+gap+$
                   strcompress(string(long(redcat[fitgal[i],0])),/remove_all)+'.fits'
     endif

; Get the distortion estimators:
     d1=sqrt(total(fp[n_par-4:n_par-3]^2))*fp[3]
     d3=sqrt(total(fp[n_par-2:n_par-1]^2))*fp[3]

; Get the flexion error sizes
     sig_psi1=sqrt(total(pe[n_par-4:n_par-3]^2))
     sig_psi3=sqrt(total(pe[n_par-2:n_par-1]^2))

; Store it all
     startpars[i,*]=[redcat[fitgal[i],0],sp]
     parfits[i,*] = [transpose(redcat[fitgal[i],0:2]),fp]
     parerrors[i,*]=[redcat[fitgal[i],0],pe]
     meritfigs[i,*]=[redcat[fitgal[i],0],chisq,dof,chisq/dof,$
                     d1,d3,sig_psi1,sig_psi3,$
                     rmag[fitgal[i]],color[fitgal[i]],$
                     niter,pegflag,failflag,badflag]
     corrmats[i,*] =[redcat[fitgal[i],0],convert_sym_matrix(corr)]


; Clear some stuff out
     delvarx,dimg
     delvarx,win

  endfor
  forloop_status,0,0,countwnum,/delete
  disp_scaled_image,0,imgwnum,/delete
  
; Save the output to file

  save_data,startpars,startpars_file,comment=startpars_com
  save_data, parfits , parfits_file ,comment=parfits_com
  save_data,parerrors,parerrors_file,comment=parerrors_com
  save_data,meritfigs,meritfigs_file,comment=meritfigs_com
  save_data,corrmats, corrmats_file, comment=corrmats_com

  aim_final_status,meritfigs,10,11,12,maxiter=maxiter

  print,'Started:  '+starttime
  print,'Finished: '+systime()
;---------------------------------------------------------------------

end
