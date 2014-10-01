pro rp_fit_catalog, red_catfile, blue_catfile, imgfile, $
                    bkg_offset=bkg_offset,$
                    tag=tag, nsave=nsave,$
                    redname=redname, bluename=bluename,$
                    testlimit=testlimit,$
                    sersic=sersic

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

; Read in the data image
full_image=mrdfits(imgfile)

if keyword_set(sersic) then n_par=13 else n_par=12

; Allow for a global background offset.
if not keyword_set(bkg_offset) then bkg_offset=median(full_image[where(full_image gt 0d)])

rmag=redcat[*,8]
ralpha=redcat[*,4]
bmag=bluecat[rmag[*,0],8]
color=bmag-rmag
bkg=redcat[*,7]+bkg_offset
n_total=n_elements(color)


; ---------- SIZE/COLOR SELECTION ----------

; We'll select galaxies to fit by size and magnitude.

fitgal=interactive_select(ralpha,rmag,aname='ALPHA',bname='F775W',$
                          prompt='Select on size/magnitude',$
                          n_sel=n_fit,not_sel=notfit,n_not_sel=n_notfit)



; ---------- PSF MODEL ESTIMATION ----------

; Do something more sophisticated later.  For now, the HST PSF is
; approximately 0.09", or 1.8 pixels in FWHM.  We'll use this
; for all the galaxies, though later we'll do something more
; individualized. 

psf_model=mk_gaussian_psf(1.8d)

; ---------- IMAGE FITTING ----------

; Make the data saving arrays.
startpars=dblarr(n_fit,n_par+1) ; Index and each starting parameter.
parfits = dblarr(n_fit,n_par+3) ; Index, object location, and each of the 
                                ; parameter fit values.
parerrors=dblarr(n_fit,n_par+1) ; Index and each parameter error.
meritfigs=dblarr(n_fit,7)       ; Index and the following figures of merit:
                                ;    chi^2=sum( (D-M)^2/M )
                                ;    DoF (degrees of freedom)
                                ;    chi^2/DoF
                                ;    N_iter (number of iterations to
                                ;            the returned fit)
                                ;    PegFlag (A binary representation
                                ;             of the pegged
                                ;             parameters, if any)
                                ;    FailFlag (binary fail marker)

; Set up file names and their save comments
startpars_file=tag+'_start.dat'
parfits_file = tag+'_fit.dat'
parerrors_file=tag+'_err.dat'
meritfigs_file=tag+'_fom.dat'

startpars_com='Starting Parameter Values: '+$
              'N, logN0, Xc, Yc, alpha, E+, Ex, '
if keyword_set(sersic) then startpars_com+='n, '
startpars_com+='g1, g2, G11, G12, G31, G32'

parfits_com = 'Fit Parameter Values: '+$
              'N, X_field, Y_field, logN0, Xc, Yc, alpha, E+, Ex, '
if keyword_set(sersic) then parfits_com+='n, '
parfits_com+='g1, g2, G11, G12, G31, G32'

parerrors_com='Parameter Error Estimates: '+$
              'N, dlogN0, dXc, dYc, dalpha, dE+, dEx, '
if keyword_set(sersic) then parerrors_com+='dn, '
parerrors_com+='dg1, dg2, dG11, dG12, dG31, dG32'

meritfigs_com='Figures of Merit: '+$
              'N, chisq, DoF, chisq/DoF, N_iter, PegFlag, FailFlag'

; Some general fitting parameters
cutrad_scale=3d                 ; Cut with a window with radius 
                                ; 3*alpha_init (from SExtractor)

min_cutrad=20d                  ; Minimum radius for fitting
max_cutrad=200d                 ; Maximum radius for fitting

maxindex=max(redcat[fitgal,0])  ; Maximum index in the fit catalog

if not keyword_set(nsave) then nsave=1

; Set up the shear lens model
center=[3130d,2865d]            ; Center of the lens in the SWARPed image
z_lens=0.183d
theta_s=0.338d/cosmo_dist(0.183) ; Scale radius (in radians) Coe et al. 2010
theta_s*=(180d/!dpi)*3600d/0.05d
c200=7.6d                        ; Also from Coe et al. 2010
a1689_lensmodel=nfw_lensmodel(redcat[*,1:2],$
                              ctr=center,theta_s=theta_s,z=z_lens,c200=c200)
shearmodel=a1689_lensmodel[*,0:1]



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
   cutrad=max(cutrad,min_cutrad)
   cutrad=min(cutrad,max_cutrad)

   win=mk_window(2d*cutrad)

; Build the PSF from a model
; For now, that model is just a uniform Gaussian
   psf=psf_model                ;Eventually this will be a function

; Cut and display the stamp
   dimg=cut_stamp(full_image,redcat[fitgal[i],1:2],win)
   disp_scaled_image,dimg-bkg[fitgal[i]]*win,imgwnum,title='Current Image',$
                     label='Object #'+$
                     strcompress(string(long(redcat[fitgal[i],0])),/remove_all)

; Clear out the start parameters and fit the image. For now
; I'll use the start parameter determination that I wrote.
   delvarx,sp
; Decide about whether to use the shear model
   sh=dblarr(2)
   if sqrt(total(shearmodel[fitgal[i],*]^2)) lt 1d then sh[*]=shearmodel[fitgal[i],*]
   fp=rp_fit_image(dimg,bkg[fitgal[i]],win,psf,$
                   startpars=sp,fixedpars=[6,7],$
                   fixedvals=sh,$
                   parerrors=pe,chisq=chisq,dof=dof,$
                   npegged=npegged,pegged=pegged,$
                   niter=niter,fail=failflag, seed=seed)
                   
; Save the data, fit and residual images if i is a multiple of nsave.
   if (i mod nsave) eq 0 then begin
      fimg=rp_mk_image(fp,redcat[fitgal[i],7],win,psf)
      rimg=fimg-dimg
      maxct=floor(alog10(maxindex))+1
      curct=floor(alog10(redcat[fitgal[i],0]))+1
      gap=''
      for j=curct,maxct-1 do gap+='0'
   
      mwrfits,dimg,tag+'_dimg'+gap+$
              strcompress(string(long(redcat[fitgal[i],0])),/remove_all)+'.fits'
      mwrfits,fimg,tag+'_fimg'+gap+$
              strcompress(string(long(redcat[fitgal[i],0])),/remove_all)+'.fits'
      mwrfits,rimg,tag+'_rimg'+gap+$
              strcompress(string(long(redcat[fitgal[i],0])),/remove_all)+'.fits'
   endif

; Create PegFlag and store the output.
   pegflag=long(total(2^pegged))

   startpars[i,*]=[redcat[fitgal[i],0],sp]
   parfits[i,*] = [transpose(redcat[fitgal[i],0:2]),fp]
   parerrors[i,*]=[redcat[fitgal[i],0],pe]
   meritfigs[i,*]=[redcat[fitgal[i],0],chisq,dof,chisq/dof,$
                   niter,pegflag,failflag]

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

print,'Done!!'
print,'N_Failed = '+strcompress(string(n_elements(where(meritfigs[*,6] gt 0))),/remove_all)
print,'N_Pegged = '+strcompress(string(n_elements(where(meritfigs[*,5] gt 0))),/remove_all)
print,'Started at:......'+starttime
print,'Finished at:.....'+systime()

;---------------------------------------------------------------------

end
