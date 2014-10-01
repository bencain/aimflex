FUNCTION RP_FIT_IMAGE, DATAIMG, BKG, WINDOW, PSF, $
                       IMG_ERROR=IMG_ERROR,$
                                ; If set, then this is the standard
                                ; deviation of each pixel.  Otherwise
                                ; the sqrt of the model value will be
                                ; used. 
                       SERSIC=SERSIC, $       
                                ; Whether to use the Sersic profile
                       MOFFAT=MOFFAT, $       
                                ; Whether to use the Moffat profile
                       PSEUDOGAUSS=PSEUDOGAUSS, $
                                ; Whether to use a pseudo-Gaussian
                                ; profile.
                       DOUBLEGAUSS=DOUBLEGAUSS,$
                                ; Whether to use a double-Gaussian
                                ; profile. 
                       STARTPARS=STARTPARS, $ 
                                ; Initial parameter guess.
                                ; Data-driven parameters will be
                                ; selected otherwise
                       FIXEDPARS=FIXEDPARS, FIXEDVALS=FIXEDVALS, $
                                ; Which parameters (by index) to fix
                                ; during fitting and the values to fix
                                ; them at. Both or neither should be
                                ; set, and both arrays should be of
                                ; the same length.  E.g., to set the
                                ; shears g1 = 0.2 and g2 = -0.1, set
                                ; FIXEDPARS=[npar-5,npar-4],
                                ; FIXEDVALS=[0.2,-0.1].  If
                                ; FIXEDPARS=-1 then there will be no
                                ; fixed parameters.
                                ;

;;;;;;;;;;;;;;;;;;;;;;;;;;;; Figures of merit for the fit;;;;;;;;;;; 
                       PARERRORS=PARERRORS, $ 
                                ; Array of parameter errors, with
                                ; value -1 for any parameters which
                                ; returned a NaN on the covmat diagonal.
                       CHISQ=CHISQ, $         
                                ; Chi^2= (dataimg-modelimg)^2/modelimg
                                ; value for the returned fit.  
                       DOF=DOF, $             
                                ; The number of fit degrees of freedom.
                       COVMAT=COVMAT, $       
                                ; The covariance matrix at the fit
                                ; parameters. The diagonal elements
                                ; provide PARERRORS. This matrix will
                                ; be NFIT x NFIT, and include only
                                ; those parameters which are fit.
                       CORRMAT=CORRMAT, $     
                                ; The rescaled parameter correlation
                                ; matrix, determined with
                                ; CORRMAT[i,j] =
                                ;    COVMAT[i,j]/(PARERRORS[i]*PARERRORS[j])
                                ; This matrix will be NFIT x NFIT, as
                                ; with COVMAT
                       DETCORR=DETCORR, $     
                                ; The determinant of corrmat.  The
                                ; further from unity, the more
                                ; correlations in the parameters. 0
                                ; will be returned if the fit fails.
                       SIGSQPRODUCT=SIGSQPRODUCT,$
                                ; The product of the diagonal elements
                                ; of the covariance matrix.  This is
                                ; an estimate of the error volume.
                       NPEGGED=NPEGGED,$
                                ; The number of parameters which are
                                ; pegged to one of the limits of the
                                ; parameter ranges.
                       PEGGED=PEGGED,$
                                ; An array of the parameters which are
                                ; pegged to one of the limits.  This
                                ; will return -1 if none are pegged.
                       PEGFLAG=PEGFLAG,$
                                ; long(total(2^pegged))
                       NITER=NITER, $         
                                ; The number of iterations required to
                                ; reach the fit.
                       FAILFLAG=FAILFLAG, $ 
                                ; As with PEGFLAG, an indicator of
                                ; which parameters had bad errors.

;;;;;;;;;;;;;;;;;;;;;;;;;;;; Other fitting settings;;;;;;;;;;; 
                       FIT_INDEX=FIT_INDEX,$
                                ; Set this to output a numeric
                                ; indicator of the fit if it fails.
                       VERBOSE=VERBOSE,$
                                ; Set this to output data on the
                                ; fitting. Not usually helpful.
                       MAXITER=MAXITER,$
                                ; The maximum iterations before the
                                ; fit just exits and returns the
                                ; current parameter set.
                       SEED=SEED
                                ; The RNG seed.  Be sure to use if
                                ; fitting more than one object.

; This is a wrapper to use MPFIT to fit for the best model parameters
; for an image.  We estimate good starting parameters, then fit the
; model to the data.  

if not keyword_set(img_error) then img_error=0d
if not keyword_set(maxiter) then maxiter=500L
if (keyword_set(sersic) or keyword_set(moffat)) then npar=13 else $
   if (keyword_set(pseudogauss) or keyword_set(doublegauss)) then npar=14 else npar=12
niter=-1

label=''
if keyword_set(fit_index) then $
   label+=strcompress(string(long(fit_index[0])),/remove_all)+' - '


dims=size(dataimg,/dimensions)

; Set parameter limits
parinf=replicate({value:0d, fixed:0, limited:[1,1], limits:[-1d,1d],$
                  parname:'',mpmaxstep:0},npar)

plimits=dblarr(npar,2)

plimits[0,*]=[-10d,10d]               ; This is a BIG range.
;plimits[1,*]=[-0.75d,0.75d]*min(dims) ; Keep Xc,Yc nearish to the center.
;plimits[2,*]=[-0.75d,0.75d]*min(dims) ;

plimits[1,*]=[-2d,2d]*min(dims) ; Wide range on Xc and Yc
plimits[2,*]=[-2d,2d]*min(dims) ;

plimits[3,*]=[1d,$                    ;
              1.5d*double(min(dims))]      ; Alpha between 1 and 3r_stamp


;plimits[4,*]=[-1d,1d]         ; We'll hope E doesn't approach unity too often
;plimits[5,*]=[-1d,1d]
; 11 May 2010 - the fitting runs to unity too often.  The problem is
;               that it goes to |E|=1-eps too often, not that the
;               parameters run over the magnitude limit.  So
;               we'll try to pare down the space.
e_lim=0.999d
plimits[4,*]=[-e_lim,e_lim]
plimits[5,*]=[-e_lim,e_lim]

g_limit=3d                      
                                ; Shear limiting.  This is a
                                ; large limit since we'll
                                ; usually fix the shear.
plimits[npar-6,*]=[-1d,1d]*g_limit 
plimits[npar-5,*]=[-1d,1d]*g_limit

psiN_limit=0.5d 
                                ; This is definitely large.
                                ; We'll likely toss any fits
                                ; with psiN anywhere near this.  
plimits[npar-4,*]=[-1d,1d]*psiN_limit
plimits[npar-3,*]=[-1d,1d]*psiN_limit
plimits[npar-2,*]=[-1d,1d]*psiN_limit
plimits[npar-1,*]=[-1d,1d]*psiN_limit

; Set the parameter names
parinf[0:5].parname=['logN0','X_c','Y_c','Alpha','E_+','E_x']
parinf[npar-6:npar-1].parname=['g_1','g_2','psi1_1','psi1_2','psi3_1','psi3_2']
pnames=parinf[*].parname

if keyword_set(sersic) then begin
; profile propto exp(-(r/alpha)^(1/n))
   plimits[6,*]=[0.2d,10d]
   parinf[6].parname='n'
endif else if keyword_set(pseudogauss) then begin
; profile propto 1/(1 + x^2 + k1*x^4/2 + k2*x^6/6)
   klim=10d
   plimits[6,*]=[0d,klim]
   parinf[6].parname='k1'
   plimits[7,*]=[0d,klim]
   parinf[7].parname='k2'
endif else if keyword_set(doublegauss) then begin
; profile is the sum of two Gaussians
   plimits[6,*]=[0d,8d]
   parinf[6].parname='logN0_2'
   plimits[7,*]=[1d-2,0.25d*double(min(dims))]
   parinf[7].parname='alpha_2'
endif else if keyword_set(moffat) then begin
; profile propto 1/(1 + x^2)^b
   plimits[6,*]=[1d,10d]
   parinf[6].parname='b'
endif

; Centers and offset ranges
pmids=(plimits[*,1]+plimits[*,0])/2d
pranges=(plimits[*,1]-plimits[*,0])/2d
for i=0,npar-1 do pranges[i]=max([pranges[i],1d-10])

; Limit the stepsize in parameter space
parinf[*].mpmaxstep=1d-3

; Define the fit and fixed parameter sets
nfit=npar
nfixed=0
if keyword_set(fixedpars) then begin
   fitpars=set_difference(lindgen(npar),fixedpars,count=nfit)
   nfixed=npar-nfit
endif else fitpars=lindgen(npar)

; If verbose is set, put out some information.
if keyword_set(verbose) then begin
   print,'****'
   for i=0,nfit-1 do $
      print,label+'Fitting for ',parinf[fitpars[i]].parname
   for i=0,nfixed-1 do $
      print,label+'Fixing '+parinf[fixedpars[i]].parname+' at '+$
            strcompress(string(fixedvals[i]),/remove_all)
   print,'****'
endif

; Set the fixed values offset by pmids and scaled by pranges
if nfixed gt 0 then begin
   parinf[fixedpars].fixed=1
   parinf[fixedpars].value=(fixedvals-pmids[fixedpars])/pranges[fixedpars]
endif


; Inputs to the fitting
fname='rp_image_deviate'
fargs={dataimg:dataimg, bkg:bkg, window:window, psf:psf, $
       sersic:keyword_set(sersic), pseudogauss:keyword_set(pseudogauss),$
       doublegauss:keyword_set(doublegauss),moffat:keyword_set(moffat),$
       pranges:pranges, pmids:pmids, fixed_error:img_error}

; Put in the starting parameters
if not keyword_set(startpars) then begin
   startpars=rp_start_pars(dataimg,bkg,nsigma=0d,seed=seed,$
                           sersic=keyword_set(sersic),$
                           pseudogauss=keyword_set(pseudogauss),$
                           doublegauss=keyword_set(doublegauss),$
                           moffat=keyword_set(moffat))
endif
if nfixed gt 0 then startpars[fixedpars]=fixedvals

; Push all parameters to a range of -1 to 1
parinf[fitpars].value=(startpars[fitpars]-pmids[fitpars])/pranges[fitpars]

pars=mpfit(fname, functargs=fargs, covar=covar, parinfo=parinf, $
           bestnorm=chisq, dof=dof, quiet=(not keyword_set(verbose)), $
           niter=niter, maxiter=maxiter,npegged=npegged)

; Check if it iterated out instead of converging
if niter ge maxiter then print,label+'RP_FIT_IMAGE.PRO: Exceeded maximum number of iterations '+$
                               strcompress('('+string(maxiter)+').',/remove_all)

; Find out which, if any parameters are pegged
if npegged gt 0 then begin
   pegged=where((pars eq 1d) or (pars eq -1d),count)
   p_hi=where(pars eq 1d,nphi)
   p_lo=where(pars eq -1d,nplo)
   pegtag=strarr(npar)
   if nphi gt 0 then pegtag[p_hi]=' (hi)'
   if nplo gt 0 then pegtag[p_lo]=' (lo)'

   print,label+'RP_FIT_IMAGE.PRO: Pegged parameters:'
   pegmsg='    '
   for j=0,count-2 do pegmsg+=parinf[pegged[j]].parname+pegtag[pegged[j]]+', '
   pegmsg+=parinf[pegged[count-1]].parname+pegtag[pegged[count-1]]
   print,pegmsg
endif else pegged=-1
pegflag=long(total(2^pegged))

; Undo the strict unity-magnitude transformation
pars=pars*pranges+pmids
   
; Get the errors from the covariance matrix, with a rescaling from the
; limited range transform
parerrors=dblarr(npar)-1d
for i=0,nfit-1 do begin
   if n_elements(covar) eq npar^2 then begin
      parerrors[fitpars[i]]=sqrt(covar[fitpars[i],fitpars[i]])*pranges[fitpars[i]] 
   endif
endfor

; Check to see if we failed somewhere.
fail=0L
baderrs=set_intersection(fitpars,where(parerrors le 0d),count=nbad)
if nbad gt 0 then begin
   fail=1L
   print,label+'RP_FIT_IMAGE.PRO: Zero error value for:'
   badmsg='    '
   for j=0,nbad-2 do badmsg+=parinf[baderrs[j]].parname+', '
   badmsg+=parinf[baderrs[nbad-1]].parname
   print,badmsg
endif
failflag=long(total(2^baderrs))


; Make and rescale the covariance and correlation matrices. These will
; be NFIT x NFIT matrices.
covmat= dblarr(nfit,nfit)
corrmat=dblarr(nfit,nfit)

for i=0,nfit-1 do for j=0,nfit-1 do begin
   if n_elements(covar) eq npar^2 then begin
      covmat[i,j]=covar[fitpars[i],fitpars[j]]*$
                  pranges[fitpars[i]]*pranges[fitpars[j]]
      corrmat[i,j]=covar[fitpars[i],fitpars[j]]/$
                   (sqrt(covar[fitpars[i],fitpars[i]]*$
                         covar[fitpars[j],fitpars[j]]) + 1d-10)
      ; The 1d-10 is there to deal with possible zero valued errors.
   endif
endfor

; Calculate DETCORR and SIGSQPRODUCT
if nbad lt 1 then begin
   detcorr=la_determ(corrmat,/check,/double,zero=1d-20)
   sigsqproduct=product(parerrors^2) 
endif else begin
   sigsqproduct=-1
   detcorr=0
endelse
      
return, pars

end
