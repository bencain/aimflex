function fit_sis,position,flexions

parinf=replicate({value:0d, fixed:0, limited:[1,1], limits:[-1d,1d],$
                  parname:'',mpmaxstep:0},3)

pmins=[5d1,$                    ;Theta_E
       1d,$                     ;r/theta_e
       -!dpi]                   ;phi

pmaxs=[3d3,$                    ;Theta_E
       3d,$                     ;r/theta_e
       !dpi]                    ;phi

pmids = (pmaxs+pmins)/2d
pranges=(pmaxs-pmins)/2d

parinf[*].parname=['theta_E','r/theta_e','phi']

sp=[0d,0d,0d]

parinf[*].value=sp

fname='sis_deviate'
fargs={flexions:flexions, pmids:pmids, pranges:pranges}

fitpars=mpfit(fname, functargs=fargs, covar=covar, parinfo=parinf, $
              bestnorm=chisq, dof=dof, /quiet,$
              niter=niter, maxiter=100,npegged=npegged)

return,fitpars*pranges+pmids

end
