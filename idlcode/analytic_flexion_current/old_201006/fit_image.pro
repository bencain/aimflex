function fit_image, dataimg, bkg, window, psf, startpars=startpars, $
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
if not keyword_set(ntry) then ntry=1
if not keyword_set(maxiter) then maxiter=250L
if not keyword_set(flexionscale) then flexionscale=1d3
if keyword_set(sersic) then npar=13 else npar=12
niter=-1

dims=size(dataimg,/dimensions)

; Inputs to the fitting
fname='image_deviate'
fargs={dataimg:dataimg, bkg:bkg, window:window, psf:psf, $
       flexionscale:flexionscale,sersic:keyword_set(sersic)}

; Set parameter limits
parinf=replicate({value:0d, fixed:0, limited:[1,1], limits:[-1d,1d],$
                  parname:'',mpmaxstep:0},npar)

parinf[0].limits=[0d,10d]
parinf[1:2].limits=[-0.25d,0.25d]*min(dims) ; Keep Xc,Yc near the center.
parinf[3].limits=[1d,0.5d*double(min(dims))] ; Alpha between 1 and r_stamp
parinf[4:5].limits=[-1d,1d] ; We'll hope E doesn't approach unity too often

g_limit=10d
parinf[npar-6:npar-5].limits=[-1d,1d]*g_limit ; Shear less than unity, if free

Gn_limit=0.5d ; This is definitely large.  We'll likely toss any fits with Gn
              ; anywhere near this

parinf[npar-4:npar-1].limits=[-1d,1d]*Gn_limit*flexionscale

if keyword_set(sersic) then begin
   parinf[6].limited=[1,1]
   parinf[6].limits=[0.25d,5d]
   parinf[6].parname='n'
endif

; Limit the stepsize in parameter space
parinf[*].mpmaxstep=1d-3

parinf[0:5].parname=['I0','X_c','Y_c','Alpha','e_+','e_x']
parinf[npar-6:npar-1].parname=['g_1','g_2','G1_1','G1_2','G3_1','G3_2']

pnames=parinf[*].parname

if keyword_set(fixedpars) then begin
   fitpars=set_difference(lindgen(npar),fixedpars)
endif else fitpars=lindgen(npar)
;;;
; fixedpars is a keyword to pass in an array of parameters to remain
; fixed at their starting values through the fitting.
if keyword_set(verbose) then print,'****'

for i=0,n_elements(fitpars)-1 do begin
   if keyword_set(verbose) then $
      print,'Fitting for ',parinf[fitpars[i]].parname
endfor

for i=0,n_elements(fixedpars)-1 do begin
   parinf[fixedpars[i]].fixed=1
   parinf[fixedpars[i]].value=fixedvals[i]
   if keyword_set(verbose) then $
      print,'Fixing '+parinf[fixedpars[i]].parname+' at '+$
            strcompress(string(parinf[fixedpars[i]].value),/remove_all)
endfor

if keyword_set(verbose) then print,'****'
;;;

; Try NTRY times to get a fit with non-zero error from the covariance matrix

startps=dblarr(npar)
fail=0L

for k=0,ntry-1 do begin
   
   if not keyword_set(startpars) then begin
      startps=start_pars(dataimg,bkg,nsigma=0.5d,seed=seed,$
                         sersic=keyword_set(sersic))
   endif else startps=startpars

   parinf[fitpars].value=startps[fitpars]



   pars=mpfit(fname, functargs=fargs, covar=covmat, parinfo=parinf,$
              bestnorm=chisq, dof=dof, /quiet, niter=niter,maxiter=maxiter,$
              ftol=tol,gtol=tol,xtol=tol)

;   pars*=scls
;   pars+=mids ; Undo the strict unity transformation

   parerrors=dblarr(npar)
   for i=0,npar-1 do begin

      if n_elements(covmat) eq npar^2 then sig=sqrt(covmat[i,i]) else sig=-1d

      parerrors[i]=sig
   endfor

   baderrs=set_intersection(fitpars,where(parerrors le 0d))
   if baderrs[0] lt 0 then nbad=0 else nbad=n_elements(baderrs)
   if nbad gt 0 then begin
      ;
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

corrmat=covmat

if n_elements(covmat) eq npar^2 then begin
   for i=0,npar-1 do begin
      for j=0,npar-1 do begin
         corrmat[i,j]/=(parerrors[i]*parerrors[j])
;         covmat[i,j]*=scls[i]*scls[j]
      endfor
   endfor
endif

         
return, pars

end
