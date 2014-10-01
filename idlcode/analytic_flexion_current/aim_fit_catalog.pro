pro aim_fit_catalog, imagefile, errorfile, bkgfile, catfile,$
                     tag=tag,testlimit=testlimit,$
                     maxiter=maxiter, nsave=nsave, min_cutrad=min_cutrad,$
                     max_cutrad=max_cutrad, cutrad_scale=cutrad_scale,$
                     psf_fwhm=psf_fwhm, hst_psf=hst_psf,$
                     pseudogauss=pseudogauss, fixk1=fixk1, fixk2=fixk2, $
                     sersic=sersic, moffat=moffat,$
;                     nfw_scale=nfw_scale, nfw_ctr=nfw_ctr,$
;                     detJacmin=detJacmin, detJacmax=detJacmax,$
                     shearmodel_pars=shearmodel_pars,$
                     freeshear=freeshear, fit_nolens=fit_nolens,$
                     fixpsi3=fixpsi3


; This function takes in the catalog file output of the
; SEX_CAT_CONVERT script, the full fits image, and an output filename
; for the fit catalog.  The fitting procedure is as follows:
;
;  1. Color/Magnitude/Size discrimination
;       - Fit galaxies, member galaxies, stars, bad objects
;  2. PSF estimation
;  3. Image Fitting

; Read in the catalogs.  The entries are:
;       N, X_field, Y_field, logS0, Alpha, eplus, ecross, Bkg, Mag,
;       MagErr, Class_star  
  cat=read_data(catfile,comment=catcom,/quiet)
; Read in the data image and the error image
  ext=0
  full_data_img=mrdfits(imagefile,ext)
  if n_elements(full_data_img) le 1 then begin
     ext++
     full_data_img=mrdfits(imagefile,ext)
  endif
  ext=0
  full_err_img=mrdfits(errorfile)
  if n_elements(full_err_img) le 1 then begin
     ext++
     full_err_img=mrdfits(errorfile,ext)
  endif
  ext=0
  full_bkg_img=mrdfits(bkgfile)
  if n_elements(full_bkg_img) le 1 then begin
     ext++
     full_bkg_img=mrdfits(bkgfile,ext)
  endif

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
  endif
  if keyword_set(fixpsi3) then begin
     fitps=set_difference(fitps,[n_par-2,n_par-1],count=n_fitps)
     fixps=set_union(fixps,[n_par-2,n_par-1],count=n_fixps)
  endif
 
  

  if n_fixps gt 0 then fixedvals=dblarr(n_fixps)

; Pull out the sizes and backgrounds.
  alpha=cat[*,4]
  Emag=sqrt(total(cat[*,5:6]^2,2))
  eps=sqrt((1d - Emag)/(1d + Emag))
  size=alpha/sqrt(eps)
  bkg=cat[*,7]
  n_fit=n_elements(alpha)


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
  meritfigs=dblarr(n_fit,12)       ; Index and the following figures of merit:
                                ;    chi^2=sum( (D-M)^2/M )
                                ;    DoF (degrees of freedom)
                                ;    chi^2/DoF
                                ;    D1=alpha*sqrt(psi11^2 + psi12^2)
                                ;    D3=alpha*sqrt(psi31^2 + psi32^2)
                                ;    sig_psi1 (error in psi11-psi12
                                ;            space)
                                ;    sig_psi3 (error in psi31-psi32
                                ;            space)
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

  fitparlist='logS0, Xc, Yc, alpha, Eplus, Ecross, '
  if keyword_set(sersic) then fitparlist+='n, '
  if keyword_set(pseudogauss) then fitparlist+='k1, k2, '
  if keyword_set(moffat) then fitparlist+='b, '
  
  fitparlist+='g_1, g_2, psi1_1, psi1_2, psi3_1, psi3_2'
  
  startpars_com='Starting parameter values: N, '+fitparlist

  parfits_com = 'Fit parameter values: N, X_field, Y_field, '+fitparlist

  parerrors_com='Fit parameter errors: N, '+fitparlist

  meritfigs_com='Figures of Merit: '+$
                'N, chisq, DoF, chisq/DoF, D1, D3, sig_psi1, sig_psi3, '+$
                'N_iter, PegFlag, FailFlag, BadFlag'

  corrmats_com='Correlation Matrices: '+$
               'N, [Converted 1D correlation matrices]'

; Some general fitting parameters
  if not keyword_set(cutrad_scale) then cutrad_scale=2d ; Window cut radius
  if not keyword_set(min_cutrad) then min_cutrad=15d    ; Minimum radius for fitting
  if not keyword_set(max_cutrad) then max_cutrad=200d   ; Maximum radius for fitting

  maxindex=max(cat[*,0]) ; Maximum index in the fit catalog

  if not keyword_set(nsave) then nsave=-1
  if not keyword_set(maxiter) then maxiter=1000

; Set up the shear lens model.  The shear lens model will depend on
; two inputs: positions and shearmodel_pars.  Everything else will be
; external. The return value will be an ngal x 2 shear array.
  if not keyword_set(shearmodel_pars) then input_shear=cat[*,1:2]*0d else $
     input_shear=aim_input_shear(cat[*,1:2],shearmodel_pars)

  if keyword_set(fit_nolens) then input_shear*=0d


; -------- Cut and fit each image ------------

  if keyword_set(testlimit) then n_fit=10
  starttime=systime()
  print,'Fitting Started: '+starttime

  clear_windows

  for i=0,n_fit-1 do begin
     if (i mod 10) eq 0 then $
        forloop_status,i,n_fit,countwnum,label='-=< Fitting the selected galaxies >=-'
   
; Create the windows for cutting the data and the errors
     se_alpha=cat[i,4]
     se_Emag= sqrt(total(cat[i,5:6]^2))
     se_eps = sqrt((1d - se_Emag)/(1d + se_Emag))
     
     cutrad=cutrad_scale*se_alpha/sqrt(se_eps)
     cutrad=max([cutrad,min_cutrad])
     cutrad=min([cutrad,max_cutrad])

     win=mk_window(2d*cutrad)
     ewin=win
     ewin[*]=1d

; Build the PSF from a model
; For now, that model is just a uniform Gaussian
     psf=psf_model              ;Eventually this will be a function

; Cut and display the stamp
     dimg=cut_stamp(full_data_img,cat[i,1:2],win)
     eimg=cut_stamp(full_err_img,cat[i,1:2],ewin)
     bimg=cut_stamp(full_bkg_img,cat[i,1:2],win) ; bkg image
     disp_scaled_image,alog(dimg-bkg[i]*win + 1d),$
                       imgwnum,title='Current Image',label='Object #'+$
                       strcompress(string(long(cat[i,0])),/remove_all)

; Decide about whether to use the shear model or to just set shear to zero.
     if (not keyword_set(freeshear)) then begin
              fixedvals[where(fixps eq n_par-6)]=input_shear[i,0] 
              fixedvals[where(fixps eq n_par-5)]=input_shear[i,1] 
     endif


; Set up parameter fixing for the Pseudo-Gaussian model
     if (keyword_set(fixk2) and keyword_set(pseudogauss)) then begin
        fixedvals[where(fixps eq 7)]=1d
     endif
     if (keyword_set(fixk1) and keyword_set(pseudogauss)) then begin
        fixedvals[where(fixps eq 6)]=1d
     endif

     delvarx,sp

;     fp=aim_fit_image(dimg,bkg[i],win,psf,$
     fp=aim_fit_image(dimg,bimg,win,psf,$
                     img_error=eimg, maxiter=maxiter,$
                     startpars=sp,fixedpars=fixps,fixedvals=fixedvals,$
                     parerrors=pe,chisq=chisq,dof=dof, seed=seed,$
                     pegflag=pegflag, niter=niter,failflag=failflag,$
                     corrmat=corr,$
                     sersic=keyword_set(sersic), $
                     pseudogauss=keyword_set(pseudogauss),$
                     moffat=keyword_set(moffat),$
                     fit_index=cat[i,0])

     badflag=((pegflag gt 0) or (niter ge maxiter) or (failflag gt 0))

; Save the data, fit and residual images if i is a multiple of nsave.
     if ((nsave gt 0) and ((i mod nsave) eq 0)) or (badflag ne 0) then begin
;        fimg=aim_mk_image(fp,bkg[i],win,psf,$
        fimg=aim_mk_image(fp,bimg,win,psf,$
                          sersic=keyword_set(sersic),$
                          pseudogauss=keyword_set(pseudogauss),$
                          moffat=keyword_set(moffat))
        rimg=dimg-fimg
        maxct=floor(alog10(maxindex))+1
        curct=floor(alog10(cat[i,0]))+1
        gap=''
        for j=curct,maxct-1 do gap+='0'
        
        mwrfits,fimg,tag+'_fimg'+gap+strcompress(string(long(cat[i,0])),/remove_all)+'.fits'
        mwrfits,rimg,tag+'_rimg'+gap+strcompress(string(long(cat[i,0])),/remove_all)+'.fits'
        mwrfits,dimg,tag+'_dimg'+gap+strcompress(string(long(cat[i,0])),/remove_all)+'.fits'
     endif

; Get the distortion estimators:
     d1=sqrt(total(fp[n_par-4:n_par-3]^2))*sp[3]
     d3=sqrt(total(fp[n_par-2:n_par-1]^2))*sp[3]

; Get the flexion error sizes
     sig_psi1=sqrt(total(pe[n_par-4:n_par-3]^2))
     sig_psi3=sqrt(total(pe[n_par-2:n_par-1]^2))

; Store it all
     startpars[i,*]=[cat[i,0],sp]
     parfits[i,*] = [transpose(cat[i,0:2]),fp]
     parerrors[i,*]=[cat[i,0],pe]
     meritfigs[i,*]=[cat[i,0],chisq,dof,chisq/dof,$
                     d1,d3,sig_psi1,sig_psi3,$
                     niter,pegflag,failflag,badflag]
     corrmats[i,*] =[cat[i,0],convert_sym_matrix(corr)]


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

  aim_final_status,meritfigs,8,9,10,maxiter=maxiter

  print,'Started:  '+starttime
  print,'Finished: '+systime()
;---------------------------------------------------------------------

end
