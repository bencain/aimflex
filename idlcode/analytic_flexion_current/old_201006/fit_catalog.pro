pro fit_catalog, red_catfile, blue_catfile, imgfile, $
                 tag=tag, nsave=nsave, im_ext=im_ext,$

                 startfix_matrix=startfix_matrix,$; N_totalx(npar+2) matrix
                 input_pars=input_pars, psftag=psftag,$
                 npsfgrid=npsfgrid, redname=redname, bluename=bluename,$
                 flexionscale=flexionscale, testlimit=testlimit,$
                 sersic=sersic, n_dev=n_dev

; This function takes in the catalog file output of the
; SEX_CAT_CONVERT script, the full fits image, and an output filename
; for the fit catalog.  The catalog is filtered and then each image is
; fit by the FIT_IMAGE function with the results saved.  Note that
; FIXEDVALS and STARTVALS are assumed to have a value for each entry
; in the catfiles, and will be filtered using the galaxy selection.

if not keyword_set(n_dev) then n_dev=-1

redcat=read_data(red_catfile,comment=rcatcom,/quiet)
bluecat=read_data(blue_catfile,comment=bcatcom,/quiet)

bkg=redcat[*,7]

mag=redcat[*,8]
magerr=redcat[*,9]
color=bluecat[*,8]-mag
colorerr=sqrt(magerr^2+bluecat[*,9]^2)
size=redcat[*,4] ; alpha=sqrt(AB)

n_total=n_elements(mag)

; Labels for the red and blue filters used.
if not keyword_set(redname) then redname="F775"
if not keyword_set(bluename) then bluename="F475"

if not keyword_set(flexionscale) then flexionscale=1d3

;Check for a specified image extension
if not keyword_set(im_ext) then fullimg=mrdfits(imgfile) else $
   fullimg=mrdfits(imgfile,im_ext)

; For saving fit images for troubleshooting
if not keyword_set(tag) then tag='fit_catalog_out'

fitfile=tag+'_fit.dat'
errfile=tag+'_err.dat'
psffile=tag+'_psf.dat'
startfile=tag+'_start.dat'

if not keyword_set(nsave) then nsave=-1

sel=0
good=lindgen(n_total)
n_good=n_total
pts=good
n_pts=n_good

;clear
sel=get_num_resp(['Select on size/color/magnitude?','(1=yes, 0=no)'],$
                 default=1)

if keyword_set(sersic) then npar=13 else npar=12
pnames=strarr(npar)
pnames[0:5]=['I0','X_c','Y_c','Alpha','E_+','E_x']
pnames[npar-6:npar-1]=['g_1','g_2','G1_1','G1_2','G3_1','G3_2']
if keyword_set(sersic) then pnames[6]='n'

; Do some input control on STARTFIX_MATRIX and INPUT_PARS

if not keyword_set(startfix_matrix) then begin
   startfix_matrix=dblarr(n_total,npar+2)
   input_pars=dblarr(n_total,npar)
endif

sfm_size=size(startfix_matrix,/dimensions)
if (sfm_size[0] ne n_total) or (sfm_size[1] ne (npar+2)) then begin
   startfix_matrix=dblarr(n_total,npar+2)
   print,'Bad startfix_matrix!!!'
endif

for i=0,n_total-1 do begin
   startfix_matrix[i,0]=total(startfix_matrix[i,2:npar+1] eq 2)
   startfix_matrix[i,1]=total(startfix_matrix[i,2:npar+1] eq 1)
endfor

if not keyword_set(input_pars) then input_pars=dblarr(n_total,npar)+1d-6
inp_size=size(input_pars,/dimensions)
if (inp_size[0] ne n_total) or (inp_size[1] ne npar) then begin
   startfix_matrix[*]=0d
   input_pars=dblarr(n_total,npar)
   print,'No input_pars'
endif


while sel do begin
; Select galaxies/stars in color/magnitude/size space.

   good=catalog_select(mag,color,size,pstag=tag,/toeps,$
                       ;/sel_upperonly,$
                       obj=obj,n_obj=n_obj, gals=gals,n_gals=n_gals,$
                       pts=pts,n_pts=n_pts, redseq=redseq,n_redseq=n_redseq,$
                       notsel=notsel, n_notsel=n_notsel, $
                       nondet=nondet, n_nondet=n_nondet, n_good=n_good,$
                       magtag=redname, colortag=bluename+' - '+redname,$
                       limits=limits)

   limcom='Selection limits: ptsrc mag (c,w), ptsrc size (c,w)'
   n_selection=n_elements(limits)
   if n_selection gt 6 then limcom+=', rs mag (c,w), rs col (c,w)'
   if (n_selection mod 2) eq 1 then limcom+=', mag limit'
   
   save_data,transpose(limits),tag+'_selectionlimits.dat',comment=limcom

   numbers=lonarr(5)

   if n_good gt 0 then begin
      numbers[0]=n_good
      print,string(numbers[0],format='(I5)')+' Usable galaxies.'
      save_data,redcat[good,*],tag+'_usable_red.cat',comment=rcatcom
      save_data,bluecat[good,*],tag+'_usable_blue.cat',comment=bcatcom
   endif

   if n_redseq gt 0 then begin
      numbers[1]=n_redseq
      print,string(numbers[1],format='(I5)')+' Red sequence galaxies.'
      save_data,redcat[redseq,*],tag+'_redseq.cat',comment=rcatcom
   endif

   if n_pts gt 0 then begin
      numbers[2]=n_pts
      print,string(numbers[2],format='(I5)')+' Point sources.'
      save_data,redcat[pts,*],tag+'_ptsrc.cat',comment=rcatcom
   endif

   if n_nondet gt 0 then begin
      numbers[3]=n_nondet
      print,string(numbers[3],format='(I5)')+' Only detected in '+redname
      save_data,redcat[nondet,*],tag+'_badcolor.cat',comment=rcatcom
   endif

   if n_notsel gt 0 then begin
      numbers[4]=n_notsel
      print,string(numbers[4],format='(I5)')+' Galaxies below the RS.'
      save_data,redcat[notsel,*],tag+'_notselected.cat',comment=rcatcom
   endif

   print,string(n_total,format='(I5)')+' Total.'

   nsum=total(numbers)
   if nsum ne n_total then print,string(n_total-nsum,format='(I5)')+$
                                ' objects got lost!!!!'

   ok=get_num_resp(['Satisfied with size/cmd selection','(1=yes, 0=no)?'],$
                   default=1)

   sel=not ok
   
   if sel then begin
      delfiles=tag+['_badcolor.cat',$
                    '_ptsrc.cat',$
                    '_redseq.cat',$
                    '_usable.cat',$
                    '_notselected.cat',$
                    '_sizehist.eps',$
                    '_sizehist.ps',$
                    '_cmd_sel.eps',$
                    '_cmd_sel.ps']

      print, 'Deleting previous selection files...'
      for i=0,n_elements(delfiles)-1 do begin
         if file_test(delfiles[i]) then begin
            spawn,'rm -f '+delfiles[i]
            print,'      '+delfiles[i]
         endif
      endfor
   endif

   clear_windows

endwhile


if keyword_set(testlimit) then if n_good gt 10 then begin
   nlimit=10
   use=randomu(seed,nlimit,/long) mod n_good
   good=use[sort(use)]
   n_good=nlimit
endif

if n_good gt 0 then begin

;Output catalogs
   outcat=dblarr(n_good,npar+8)

   errcom='Errors: N, X, Y, DI0, DXc, DYc, DAlpha, Deplus, Decross, '+$
          'Dg1, Dg2, DG1_1, DG1_2, DG3_1, DG3_2'
   errcat=dblarr(n_good,npar+3)

   startcat=dblarr(n_good,npar+3)
   startcat[*,0:2]=redcat[good,0:2]
   
   psfcom='PSF model: X, Y, {PSF_model_pars}'

; Fit the PSF model, or use the TinyTim PSF model

   usett=get_num_resp('Use TinyTim for the PSF (HST data)? (1=yes, 0=no)',$
                      default=0)

   if usett then begin
      if not keyword_set(psftag) then psftag=tag+'_a1689psf'
      psffiles=tinytim_psf(redcat,good,psftag)
   endif else begin
      if not keyword_set(npsfgrid) then npsfgrid=10
      psfmodel=fit_psfmodel(redcat[pts,*],npsfgrid,size(fullimg,/dimensions),$
                            npar=npsfpar,/circ)
   endelse


; Fit the images

   for i=0, n_good-1 do begin
      forloop_status,i,n_good,win,label='Fitting galaxies'

; Set up the starting parameters
      startps=dblarr(npar)
      startps[0]=redcat[good[i],3]
      startps[3:5]=redcat[good[i],4:6]
      if keyword_set(sersic) then startps[6]=0.5d

      plug_in_vals=where(startfix_matrix[good[i],2:npar+1] gt 0,n_val)
      if n_val gt 0 then startps[plug_in_vals]=$
         input_pars[good[i],plug_in_vals]

      startcat[i,3:npar+2]=startps

; Pick out the stamp size, with a minimum radius
      cutrad_scale=2d
      mincutrad=15d
      maxcutrad=100d
      cutrad=double(ceil(cutrad_scale*startps[3]*$
                         ((1d +sqrt(total(startps[4:5]^2)))/$
                          (1d -sqrt(total(startps[4:5]^2))))^0.25d))
      cutrad=max([cutrad,mincutrad])
      cutrad=min([cutrad,maxcutrad])
;      clear

; Do some output
      print,'Obj # = ',i+1
      print,'Location = '+strcompress(string(redcat[good[i],1])+','+$
                                      string(redcat[good[i],2]),/remove_all)
      print,'Cutrad = ',cutrad
      print,'Start pars = '
      if keyword_set(sersic) then begin
         print,startps[0:3]
         print,startps[4:6]
         print,startps[7:8]
         print,startps[9:12]
      endif else begin
         print,startps[0:3]
         print,startps[4:7]
         print,startps[8:11]
      endelse
      

      window=mk_window(ceil((2d)*cutrad))
      stamp=cut_stamp(fullimg,redcat[good[i],1:2],window)

; Display the stamp

      disp_scaled_image,stamp-bkg[good[i]]*window,disp_wnum

; Build the PSF
      if usett then begin
         psf=mrdfits(psffiles[i])
      endif else psf=mk_psfmodel(redcat[good[i],1:2],psfmodel,npsfpar)
      
      fps=where(startfix_matrix[good[i],2:npar+1] eq 1)
      
; Fit the image, either once with native error estimates or n_dev
; times and return the average and stddev.

      if n_dev le 3 then begin
         if startfix_matrix[good[i],1] gt 0 then begin
            fitpars=$
               fit_image(stamp,bkg[good[i]],window,psf,startpars=startps,$
                         chisq=chisq,dof=dof, parerrors=errs, $
                         flexionscale=flexionscale,$
                         fixedpars=fps,$
                         fixedvals=input_pars[good[i],fps],$
                         covmat=covmat,niter=niter,$
                         sersic=keyword_set(sersic))
         endif else begin
            fitpars=$
               fit_image(stamp,bkg[good[i]],window,psf,startpars=startps,$
                         chisq=chisq,dof=dof, parerrors=errs,covmat=covmat,$
                         niter=niter,$
                         flexionscale=flexionscale,$
                         sersic=keyword_set(sersic))
         endelse
      endif else begin
         dev_fits=dblarr(n_dev,npar)
         for j=0,n_dev-1 do begin
;            forloop_status,j,n_dev,devwin
            fail=0L

            if startfix_matrix[good[i],1] gt 0 then begin
               fitpars=$
                  fit_image(stamp,bkg[good[i]],window,psf,startpars=startps,$
                            chisq=chisq,dof=dof, parerrors=errs, $
                            flexionscale=flexionscale,$
                            fixedpars=fps,fail=fail,$
                            fixedvals=input_pars[good[i],fps],$
                            covmat=covmat,niter=niter,$
                            sersic=keyword_set(sersic))
            endif else begin
               fitpars=$
                  fit_image(stamp,bkg[good[i]],window,psf,startpars=startps,$
                            chisq=chisq,dof=dof, parerrors=errs,$
                            covmat=covmat,fail=fail,$
                            niter=niter,$
                            flexionscale=flexionscale,$
                            sersic=keyword_set(sersic))
            endelse

            dev_fits[j,*]=fitpars
            if fail then break
         endfor
;         forloop_status,0,0,devwin,/delete
;         delvarx,devwin
; Average the deviates and get errors provided that we didn't fail.
         if not fail then for j=0,npar-1 do begin
            fitpars[j]=mean(dev_fits[*,j])
            errs[j]=stddev(dev_fits[*,j])
         endfor 

      endelse


      print,'Fitpars = '
      if not keyword_set(sersic) then begin
         print,fitpars[0:3]
         print,fitpars[4:7]
         print,fitpars[8:11]
      endif else begin
         print,fitpars[0:3]
         print,fitpars[4:6]
         print,fitpars[7:8]
         print,fitpars[9:12]
      endelse

      print,'Chisq, dof = ',chisq,dof
      print,'#####################################'


      parstring=''
      for d=0,npar-2 do parstring+=pnames[d]+', '
      parstring+=pnames[npar-1]

      if (nsave gt 0) and (((i+1) mod nsave) eq 0) then begin
         if good[i]+1 lt 10 then fill='000' else $
            if good[i]+1 lt 100 then fill='00' else $
               if good[i]+1 lt 1000 then fill='0' else fill=''
         fitimg=mk_image(fitpars,bkg[good[i]],window,psf)
         
         label=fill+strcompress(string(good[i]+1),/remove_all)

         mwrfits,stamp,tag+'_data_'+label+'.fits'
         mwrfits,fitimg,tag+'_fit_'+label+'.fits'
         mwrfits,stamp-fitimg,tag+'_resid_'+label+'.fits'
         save_data,covmat,tag+'_covmat_'+label+'.dat',$
                   comment='Covariance matrix: '+parstring
      endif

      outcat[i,0:2]=redcat[good[i],0:2]
      outcat[i,3:npar+2]=fitpars
      outcat[i,npar+3:npar+4]=[chisq,dof]
      outcat[i,npar+5]=double(chisq)/double(dof)
      outcat[i,npar+6]=niter
      outcat[i,npar+7]=bkg[good[i]]
      
      errcat[i,0:2]=redcat[good[i],0:2]
      errcat[i,3:npar+2]=errs

; Clear things out, just to be sure
      stamp*=0d
      window*=0d

   endfor

   catstring='N, X, Y, '

   outcom=catstring+parstring+', chisq, dof, chisq/dof, n_iter, bkg'
   startcom='Start parameters'+catstring+parstring

   forloop_status,0,0,win,/delete
   disp_scaled_image,0,disp_wnum,/delete

   if usett then begin
      save_data,outcat,fitfile,comment=outcom
      save_data,errcat,errfile,comment=errcom
   endif else begin
      save_data,psfmodel,psffile,comment=psfcom
      save_data,outcat,fitfile,comment=outcom
      save_data,errcat,errfile,comment=errcom
   endelse

   save_data,startcat,startfile,comment=startcom

endif

print,'Done at '+systime()

end
