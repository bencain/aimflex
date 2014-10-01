function mcmc_shearfield_update, LENSMODEL_NOW, WEIGHTS, $
                                 FLEXIONSCALE=FLEXIONSCALE, $
                                 ERRS=ERRS, CHISQ=CHISQ, SEED=SEED

; LENSMODEL_NOW is an Nx8 array with the current lens parameters.
; Columns 0-1 are the x and y positions.  Columns 2-3 are the current
; real and imaginary shear components, respectively.  Columns 4-5 are
; the current 1-flexion and columns 6-7 are the current 3-flexion.
; WEIGHTS is an NxN array of the pair weightings.

n=(size(lensmodel_now,/dimensions))[0]
if not keyword_set(flexionscale) then flexionscale=1d3

shears=dblarr(2*n)
shears[0:n-1]=lensmodel_now[*,2]
shears[n:2*n-1]=lensmodel_now[*,3]

chisq=shearfield_chisq(shears,lensmodel_now=lensmodel_now, weights=weights,$
                       flexionscale=flexionscale)

g_max=1d
g_min=-1d

stepsize=(g_max-g_min)/100d

n_mcmc=1d3
n_save=50d

save_steps=dblarr(n_save,2*n+1)
; We'll save the best n_save steps. It will be in the format (chisq,[shears]) 
save_steps[0,*]=[chisq,shears]
chisq_split=dblarr(n_save,2)

; Make sure that we fill the save array
for i=1,n_save-1 do begin

   test_shears=shears+stepsize*randomn(seed,2*n)
   chisq_test=shearfield_chisq(test_shears,lensmodel_now=lensmodel_now,$
                               weights=weights,flexionscale=flexionscale,$
                               chisq_diff=chisq_diff, chisq_norm=chisq_norm)
;   print,chisq_test
   accept_limit=exp(chisq-chisq_test)
   accept_query=randomu(seed,/double)

   if accept_query lt accept_limit then begin
      chisq=chisq_test
      shears=test_shears
   endif

   save_steps[i,*]=[chisq,shears]
   chisq_split[i,*]=[chisq_diff,chisq_norm]
endfor

order=sort(save_steps[*,0])
save_steps=save_steps[order,*]

for i=n_save,n_mcmc-1 do begin

   if ((i+1) mod n_save) eq 0 then forloop_status,i,n_mcmc,wnum

   test_shears=shears+stepsize*randomn(seed,2*n)
   chisq_test=shearfield_chisq(test_shears,lensmodel_now=lensmodel_now,$
                               weights=weights,flexionscale=flexionscale,$
                               chisq_diff=chisq_diff, chisq_norm=chisq_norm)

   accept_limit=min([1d,exp(chisq-chisq_test)])
   accept_query=randomu(seed,/double)

   if accept_query lt accept_limit then begin
      chisq=chisq_test
      shears=test_shears
   endif

   if chisq lt save_steps[n_save-1,0] then begin
      save_steps[n_save-1,*]=[chisq,shears]
      order=sort(save_steps[*,0])
      save_steps=save_steps[order,*]
      chisq_split[n_save-1,*]=[chisq_diff,chisq_norm]
      chisq_split=chisq_split[order,*]
   endif

;   save_steps[i mod n_save,*]=[chisq,shears]
;   chisq_split[i mod n_save,*]=[chisq_diff,chisq_norm]
endfor

forloop_status,0,0,wnum,/delete

best_shear=dblarr(n,2)
errs=dblarr(n,2)

;for i=0,n-1 do best_shear[i,0]=mean(save_steps[*,i+1])
;for i=0,n-1 do best_shear[i,1]=mean(save_steps[*,i+1+n])

chisq=min(save_steps[*,0],best)
best_shear[*,0]=save_steps[best,1:n]
best_shear[*,1]=save_steps[best,n+1:2*n]

;chisq=shearfield_chisq([best_shear[*,0],best_shear[*,1]],$
;                       lensmodel_now=lensmodel_now,$
;                       weights=weights,flexionscale=flexionscale)

for i=0,n-1 do errs[i,0]=stddev(save_steps[*,i+1])
for i=0,n-1 do errs[i,1]=stddev(save_steps[*,i+1+n])

window,xsize=900,ysize=600
!p.multi=[0,2,2]
plot,save_steps[*,0],title='Chi-squared'
plot,save_steps[*,1],title='galaxy 2 g_1'
plot,save_steps[*,1+n],title='galaxy 2 g_2'
plot,chisq_split[*,0],title='Chisq breakdown'
oplot,chisq_split[*,1],linestyle=2
!p.multi=0

return, best_shear

end
