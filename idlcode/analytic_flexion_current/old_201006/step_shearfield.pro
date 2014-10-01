function step_shearfield, lensmodel_now, weights, $
                          flexionscale=flexionscale, $
                          ntry=ntry, maxiter=maxiter, parerrors=parerrors, $
                          chisq=chisq, dof=dof, covmat=covmat, $
                          corrmat=corrmat, fail=fail, niter=niter

if not keyword_set(flexionscale) then flexionscale=1d3
n=(size(gal_pos,/dimensions))[0]

if not keyword_set(ntry) then ntry=10
if not keyword_set(maxiter) then maxiter=250L

fname='shearfield_deviate'
fargs={lensmodel_now:lensmodel_now, weights:weights, $
       flexionscale:flexionscale}

parinf=replicate({value:0d, fixed:0, limited:[1,1], limits:[-0.5d,0.5d],$
                  parname:''},2*n)

startpars=dblarr(2*n)
covmat=0d

startpars[*]=0d


new_shears=mpfit(fname,functargs=fargs,bestnorm=chisq,$
                 dof=dof,niter=niter,maxiter=maxiter)

return,new_shears

end
