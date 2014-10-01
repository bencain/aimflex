function new_fit_image, dataimg, bkg, window, psf, startpars=startpars, $
                    parerrors=parerrors, chisq=chisq, covmat=covmat, $
                    corrmat=corrmat, dof=dof, fixedpars=fixedpars, $
                    fixedvals=fixedvals, verbose=verbose,$
                    fail=fail,niter=niter, tol=tol, maxiter=maxiter,$
                    seed=seed, ntry=ntry, flexionscale=flexionscale,$
                    sersic=sersic

; This is a wrapper to use mpfit to fit for the best model parameters
; for an image.  We estimate good starting parameters, then fit the
; model to the data.  The output parameters are given in the POL
; format, and input STARTPARS is assumed to be in that format as well.

if not keyword_set(tol) then tol=1D-10
if not keyword_set(ntry) then ntry=3
if not keyword_set(maxiter) then maxiter=250L
if not keyword_set(flexionscale) then flexionscale=1d3
if keyword_set(sersic) then npar=13 else npar=12
niter=-1

dims=size(dataimg,/dimensions)

; Set parameter limits
parinf=replicate({value:0d, fixed:0, limited:[1,1], limits:[-1d,1d],$
                  parname:'',mpmaxstep:0},npar)

plimits=dblarr(npar,2)

plimits[0,*]=[0d,10d]
plimits[1,*]=[-0.25d,0.25d]*min(dims) ; Keep Xc,Yc near the center.
plimits[2,*]=[-0.25d,0.25d]*min(dims)
plimits[3,*]=[1d,0.5d*double(min(dims))] ; Alpha between 1 and r_stamp
plimits[4,*]=[-1d,1d] ; We'll hope E doesn't approach unity too often
plimits[5,*]=[-1d,1d]

g_limit=1d  ; Shear less than unity, if free
plimits[npar-6,*]=[-1d,1d]*g_limit 
plimits[npar-5,*]=[-1d,1d]*g_limit

Gn_limit=0.1d ; This is definitely large.  We'll likely toss any fits with Gn
              ; anywhere near this
plimits[npar-4,*]=[-1d,1d]*Gn_limit*flexionscale
plimits[npar-3,*]=[-1d,1d]*Gn_limit*flexionscale
plimits[npar-2,*]=[-1d,1d]*Gn_limit*flexionscale
plimits[npar-1,*]=[-1d,1d]*Gn_limit*flexionscale

if keyword_set(sersic) then begin
   plimits[6,*]=[0.2d,10d]
   parinf[6].parname='n'
endif

; Centers and offset ranges
pmids=(plimits[*,1]+plimits[*,0])/2d
pranges=(plimits[*,1]-plimits[*,0])/2d

;print,pmids
;print,pranges

; Limit the stepsize in parameter space
parinf[*].mpmaxstep=1d-6

parinf[0:5].parname=['I0','X_c','Y_c','Alpha','E_+','E_x']
parinf[npar-6:npar-1].parname=['g_1','g_2','G1_1','G1_2','G3_1','G3_2']

pnames=parinf[*].parname

; Define the fit and fixed parameter sets
nfit=npar
nfixed=0
if keyword_set(fixedpars) then begin
   fitpars=set_difference(lindgen(npar),fixedpars,count=nfit)
   nfixed=npar-nfit
endif else fitpars=lindgen(npar)

; fixedpars is a keyword to pass in an array of parameters to remain
; fixed at their starting values through the fitting.

if keyword_set(verbose) then begin
   print,'****'
   for i=0,nfit-1 do $
      print,'Fitting for ',parinf[fitpars[i]].parname
   for i=0,nfixed-1 do $
      print,'Fixing '+parinf[fixedpars[i]].parname+' at '+$
            strcompress(string(fixedvals[i]),/remove_all)
   print,'****'
endif

; Set the fixed values offset by pmids and scaled by pranges
for i=0,nfixed-1 do begin
   parinf[fixedpars[i]].fixed=1
   parinf[fixedpars[i]].value=$
      (fixedvals[i]-pmids[fixedpars[i]])/pranges[fixedpars[i]]
endfor

;;;

; Try NTRY times to get a fit with non-zero error from the covariance matrix

startps=dblarr(npar)
fail=0L

; Inputs to the fitting
fname='new_image_deviate'
fargs={dataimg:dataimg, bkg:bkg, window:window, psf:psf, $
       flexionscale:flexionscale,sersic:keyword_set(sersic),$
       pranges:pranges,pmids:pmids}


for k=0,ntry-1 do begin
   
   if not keyword_set(startpars) then begin
      startps=start_pars(dataimg,bkg,nsigma=0.5d,seed=seed,$
                         sersic=keyword_set(sersic))
   endif else startps=startpars

   ; Push all parameters to a range of -1 to 1
   parinf[fitpars].value=(startps[fitpars]-pmids[fitpars])/pranges[fitpars]

   pars=mpfit(fname, functargs=fargs, covar=covmat, parinfo=parinf,$
              bestnorm=chisq, dof=dof, /quiet, niter=niter,maxiter=maxiter,$
              ftol=tol,gtol=tol,xtol=tol)

   pars=pars*pranges+pmids ; Undo the strict unity-magnitude transformation
   
; Get the errors from the covariance matrix, with a rescaling from the
; limited range transform
   parerrors=dblarr(npar)
   for i=0,npar-1 do begin
      if n_elements(covmat) eq npar^2 then begin
         parerrors[i]=sqrt(covmat[i,i])*pranges[i] 
      endif else parerrors[i]=-1d
   endfor

; Check to see if we failed somewhere.
   baderrs=set_intersection(fitpars,where(parerrors le 0d),count=nbad)
   if nbad gt 0 then begin
      print,'FIT_IMAGE.PRO: Zero error value for:'
      badmsg='    '
      for j=0,nbad-2 do badmsg+=parinf[baderrs[j]].parname+', '
      badmsg+=parinf[baderrs[nbad-1]].parname
      print,badmsg
      print,' -- Try '+strcompress(string(k+1)+'/'+string(ntry),/remove_all)
   endif else break

   if k eq (ntry-1) then fail=1L

endfor

startpars=startps

; Make and rescale the covariance and correlation matrices
corrmat=covmat
if n_elements(covmat) eq npar^2 then begin
   for i=0,npar-1 do for j=0,npar-1 do begin
      corrmat[i,j]/=(parerrors[i]*parerrors[j])
      covmat[i,j]*=pranges[i]*pranges[j]
   endfor
endif
         
return, pars

end
