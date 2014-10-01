pro rp_fit_catalog, red_catfile, blue_catfile, imgfile, $
                    img_error=img_error, maxiter=maxiter,$
                    psf_fwhm=psf_fwhm,$
                    tag=tag, nsave=nsave,$
                    save_ximg=save_ximg, save_dimg=save_dimg,$
                    redname=redname, bluename=bluename,$
                    testlimit=testlimit, nfw_scale=nfw_scale,$
                    pseudogauss=pseudogauss, fixk1=fixk1, fixk2=fixk2, $
                    sersic=sersic, doublegauss=doublegauss, moffat=moffat,$
                    startflex=startflex, shearlim=shearlim, freeshear=freeshear

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

; Read in the data image and make sure that it's sane
  full_image=mrdfits(imgfile)
  negpix=where(full_image le 0d,complement=pospix)
  full_image[negpix]=0d

; Set up lists of the fit parameters and the fixed parameters
  if (keyword_set(sersic) or keyword_set(moffat)) then n_par=13 else $
     if (keyword_set(pseudogauss) or $
         keyword_set(doublegauss)) then n_par=14 else n_par=12
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
  fixedvals=dblarr(n_fixps)

; Pull out the magnitudes and backgrounds.  I should put zeropoints in
; at some point here.
  rmag=redcat[*,8]
  ralpha=redcat[*,4]
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

  fitgal=rp_cmr_select(color,rmag,ralpha,n_fit_obj=n_fit)


; ---------- PSF MODEL ESTIMATION ----------

; Do something more sophisticated later.  For now, the HST PSF is
; approximately 0.09", or 1.8 pixels in FWHM.  We'll use this
; for all the galaxies, though later we'll do something more
; individualized. 

  if not keyword_set(psf_fwhm) then psf_fwhm=1.8d
  psf_model=mk_gaussian_psf(psf_fwhm)

; ---------- IMAGE FITTING ----------

; Make the data saving arrays.
  startpars=dblarr(n_fit,n_par+1) ; Index and each starting parameter.
  parfits = dblarr(n_fit,n_par+3) ; Index, object location, and each of the 
                                ; parameter fit values.
  parerrors=dblarr(n_fit,n_par+1) ; Index and each parameter error.
  meritfigs=dblarr(n_fit,9)       ; Index and the following figures of merit:
                                ;    chi^2=sum( (D-M)^2/M )
                                ;    DoF (degrees of freedom)
                                ;    chi^2/DoF
                                ;    D1=alpha*sqrt(psi11^2 + psi12^2)
                                ;    D3=alpha*sqrt(psi31^2 + psi32^2)
                                ;    N_iter (number of iterations to
                                ;            the returned fit)
                                ;    PegFlag (A binary representation
                                ;             of the pegged
                                ;             parameters, if any)
                                ;    FailFlag (binary fail marker)
  n2d=n_fitps
  corrmats=dblarr(n_fit,(n2d^2-n2d)/2 + 1) 
                                ; The number of parameters -2 assumes
                                ; that the shears are fixed, and if k1
                                ; and/or k2 are fixed, subtract one
                                ; for each, since CORRMAT is returned
                                ; from RP_FIT_IMAGE as an n2d x n2d
                                ; matrix, where n2d is the number of
                                ; free parameters.

; Set up file names and their save comments
  startpars_file=tag+'_start.dat'
  parfits_file = tag+'_fit.dat'
  parerrors_file=tag+'_err.dat'
  meritfigs_file=tag+'_fom.dat'
  corrmats_file =tag+'_corr.dat'

  fitparlist='logN0, Xc, Yc, alpha, Eplus, Ecross, '
  if keyword_set(sersic) then fitparlist+='n, '
  if keyword_set(pseudogauss) then fitparlist+='k1, k2, '
  if keyword_set(doublegauss) then fitparlist+='logN0_2, alpha_2, '
  if keyword_set(moffat) then fitparlist+='b, '
  fitparlist+='g_1, g_2, psi1_1, psi1_2, psi3_1, psi3_2'

  startpars_com='Starting parameter values: N, '+fitparlist
  parfits_com = 'Fit parameter values: N, X_field, Y_field, '+fitparlist
  parerrors_com='Fit parameter errors: N, '+fitparlist

  meritfigs_com='Figures of Merit: '+$
                'N, chisq, DoF, chisq/DoF, D1, D3, N_iter, PegFlag, FailFlag'

  corrmats_com='Correlation Matrices: '+$
               'N, [Converted 1D correlation matrices]'

; Some general fitting parameters
  cutrad_scale=3d               ; Cut with a window with radius 
                                ; 3*alpha_init (from SExtractor)

  min_cutrad=4d                ; Minimum radius for fitting
  max_cutrad=200d               ; Maximum radius for fitting

                                ; Get rid of any places where galaxies
                                ; are off the edge
  fitgal=set_intersection(fitgal,where(bkg gt 0d),count=n_fit)

  maxindex=max(redcat[fitgal,0]) ; Maximum index in the fit catalog

  if not keyword_set(nsave) then nsave=-1
  if not keyword_set(maxiter) then maxiter=500

; Set up the shear lens model
  if not keyword_set(nfw_scale) then nfw_scale=1d
  center=[3130d,2865d]          ; Center of the lens in the SWARPed image

  center=[1d3,1d3] ; Just to try the dependence of the fit flexion on the lens model

  z_lens=0.183d

  r200=2.16d                   ; in h_70^-1 Mpc (Limousin et al. 2007)
  d_ang=cosmo_dist(0.183)       ; Distance to A1689 in h_70^-1 Mpc
  theta200=(r200/d_ang)*(180d/!dpi)*3600d/0.05d ; in HST pixels
  c200=7.6d                                       ; also Limousin et al. 2007
  theta_s=theta200/c200
  a1689_lensmodel=nfw_lensmodel(redcat[*,1:2],$
                                ctr=center,theta_s=theta_s,z=z_lens,c200=c200,$
                                force_scale=nfw_scale)
  if not keyword_set(shearlim) then shearlim=1d 
                                ; Limit on when to not use the shear
                                ; model and instead fix the shear to
                                ; zero. 


; -------- Cut and fit each image ------------

  if keyword_set(testlimit) then n_fit=10
  starttime=systime()
  print,'Fitting Started: '+starttime

  clear_windows

  for i=0,n_fit-1 do begin
     forloop_status,i,n_fit,countwnum,$
                    label='--=< Fitting the selected galaxies >=--'
   
; Create the window
     cutrad=cutrad_scale*redcat[fitgal[i],4]
     cutrad=max([cutrad,min_cutrad])
     cutrad=min([cutrad,max_cutrad])

     win=mk_window(2d*cutrad)

; Build the PSF from a model
; For now, that model is just a uniform Gaussian
     psf=psf_model              ;Eventually this will be a function

; Cut and display the stamp
     dimg=cut_stamp(full_image,redcat[fitgal[i],1:2],win)
     disp_scaled_image,dimg-bkg[fitgal[i]]*win,imgwnum,title='Current Image',$
                       label='Object #'+$
                       strcompress(string(long(redcat[fitgal[i],0])),/remove_all)

; Set up the starting parameters
     sp=rp_start_pars(dimg,bkg[fitgal[i]],nsigma=0d,seed=seed,$
                      sersic=keyword_set(sersic),$
                      pseudogauss=keyword_set(pseudogauss),$
                      doublegauss=keyword_set(doublegauss),$
                      moffat=keyword_set(moffat))


; Decide about whether to use the shear model or to just set shear to zero.
     if not keyword_set(freeshear) then begin
        if sqrt(total(a1689_lensmodel[fitgal[i],0:1]^2)) lt shearlim then begin
           fixedvals[n_fixps-2:n_fixps-1]=a1689_lensmodel[fitgal[i],0:1]
           sp[n_par-6:n_par-5]=fixedvals[n_fixps-2:n_fixps-1]
        endif else begin
           fixedvals[n_fixps-2:n_fixps-1]=0d
           sp[n_par-6:n_par-5]=0d
        endelse
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

     fp=rp_fit_image(dimg,bkg[fitgal[i]],win,psf,$
                     img_error=img_error, maxiter=maxiter,$
                     startpars=sp,fixedpars=fixps,fixedvals=fixedvals,$
                     parerrors=pe,chisq=chisq,dof=dof, seed=seed,$
                     pegflag=pegflag, niter=niter,failflag=failflag,$
                     corrmat=corr,$
                     sersic=keyword_set(sersic), $
                     pseudogauss=keyword_set(pseudogauss),$
                     doublegauss=keyword_set(doublegauss),$
                     moffat=keyword_set(moffat),$
                     fit_index=redcat[fitgal[i],0])
                   
; Save the data, fit and residual images if i is a multiple of nsave.
     if ((nsave gt 0) and ((i mod nsave) eq 0)) then begin
        fimg=rp_mk_image(fp,bkg[fitgal[i]],win,psf,$
                         sersic=keyword_set(sersic),$
                         doublegauss=keyword_set(doublegauss),$
                         pseudogauss=keyword_set(pseudogauss),$
                         moffat=keyword_set(moffat))
        rimg=fimg-dimg
        ximg=rimg^2/fimg
        ximg[where(fimg eq 0)]=0d
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

; Store it all
     startpars[i,*]=[redcat[fitgal[i],0],sp]
     parfits[i,*] = [transpose(redcat[fitgal[i],0:2]),fp]
     parerrors[i,*]=[redcat[fitgal[i],0],pe]
     meritfigs[i,*]=[redcat[fitgal[i],0],chisq,dof,chisq/dof,$
                     d1,d3,niter,pegflag,failflag]
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

  rp_final_status,meritfigs,6,7,8,maxiter=maxiter

  print,'Started:  '+starttime
  print,'Finished: '+systime()
;---------------------------------------------------------------------

end
