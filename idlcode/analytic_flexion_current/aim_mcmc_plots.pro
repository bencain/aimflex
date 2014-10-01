pro aim_mcmc_plots, fitfile, fomfile, trufile, errfile, stafile, tag=tag

if not keyword_set(tag) then tag='aim_mcmc_plots_' else tag+='_'

fit=read_data(fitfile,/quiet)
fom=read_data(fomfile,/quiet)
tru=read_data(trufile,/quiet)
err=read_data(errfile,/quiet)
sta=read_data(stafile,/quiet)

nfit=12
ntru=(size(tru,/dimensions))[1]-1

ok=($
   where((fom[*,3] lt 1.5d) and (fom[*,11] eq 0) and $
         (fit[*,1] gt 0d) and (fit[*,4] gt 1.8d) and $
         (fom[*,6] gt 1d-4) and (fom[*,6] lt 1d-3),nok)$
   )

print,'Nok',nok
print,'Mean/median Niter',$
      mean(fom[ok,8]),median(fom[ok,8])
print,'Mean/median chisq/DoF',$
      mean(fom[ok,3]),median(fom[ok,3])
print,'Mean/median DoF',$
      mean(fom[ok,2]),median(fom[ok,2])

print,'Mean shape errors:'
print,'logS0:  ',mean(fit[ok,1]-tru[ok,1])
print,'X_c:    ',mean(fit[ok,2]-tru[ok,2])
print,'Y_c:    ',mean(fit[ok,3]-tru[ok,3])
print,'alpha:  ',mean(fit[ok,4]-tru[ok,4])
print,'E_+:    ',mean(fit[ok,5]-tru[ok,5])
print,'E_x:    ',mean(fit[ok,6]-tru[ok,6])

print,'Shape deviation:'
print,'logS0:  ',stddev(fit[ok,1]-tru[ok,1])
print,'X_c:    ',stddev(fit[ok,2]-tru[ok,2])
print,'Y_c:    ',stddev(fit[ok,3]-tru[ok,3])
print,'alpha:  ',stddev(fit[ok,4]-tru[ok,4])
print,'E_+:    ',stddev(fit[ok,5]-tru[ok,5])
print,'E_x:    ',stddev(fit[ok,6]-tru[ok,6])


print,'Mean flexion errors:'
print,'Psi11:  ',mean(fit[ok,nfit-3]-tru[ok,ntru-3])/0.05d
print,'Psi12:  ',mean(fit[ok,nfit-2]-tru[ok,ntru-2])/0.05d
print,'Psi31:  ',mean(fit[ok,nfit-1]-tru[ok,ntru-1])/0.05d
print,'Psi32:  ',mean(fit[ok,nfit-0]-tru[ok,ntru-0])/0.05d

print,'Flexion deviation:'
print,'Psi11:  ',stddev(fit[ok,nfit-3]-tru[ok,ntru-3])/0.05d
print,'Psi12:  ',stddev(fit[ok,nfit-2]-tru[ok,ntru-2])/0.05d
print,'Psi31:  ',stddev(fit[ok,nfit-1]-tru[ok,ntru-1])/0.05d
print,'Psi32:  ',stddev(fit[ok,nfit-0]-tru[ok,ntru-0])/0.05d


vars=['logS0','X_c','Y_c','alpha','E_+','E_x','Psi_11','Psi_12','Psi_31','Psi_32']
files=tag+vars+'.ps'
xts='True '+vars
yts='Fit '+vars
;yts='Fit - True '+vars


;;;;;;;
; Fit vs true for shapes
for i=1,6 do begin
   
   setup_ps,files[i-1]
   
   plot,tru[ok,i],fit[ok,i],xtitle=xts[i-1],ytitle=yts[i-1],$
        psym=7,charsize=2
   oploterr,tru[ok,i],fit[ok,i],err[ok,i]
   oplot,[-100,100],[-100,100],color=fsc_color('red')
;   plot,tru[ok,i],fit[ok,i]-tru[ok,i],xtitle=xts[i-1],ytitle=yts[i-1],$
;        psym=7,charsize=1.5
;   oploterr,tru[ok,i],fit[ok,i]-tru[ok,i],err[ok,i]
;   oplot,[-100,100],[0,0],color=fsc_color('red')


   setup_x
;   ps2eps,files[i-1],/delete

;   plot,tru[ok,i],err[ok,i],xtitle=vars[i-1],/psym
;   wait,4


endfor

;;;;;;;;
; Fit vs true for flexions
for i=0,3 do begin

   
   setup_ps,files[i+6]
   

if i gt 1 then begin
   plot,tru[ok,i+ntru-3]/0.05d,fit[ok,i+nfit-3]/0.05d,xtitle=xts[i+6],ytitle=yts[i+6],$
        psym=7,charsize=2,xstyle=1,ystyle=1,$
        xrange=[-0.025d,0.025d],$
        yrange=[-0.05d,0.05d],$
        xticks=4,yticks=4
   oploterr,tru[ok,i+ntru-3]/0.05d,fit[ok,i+nfit-3]/0.05d,err[ok,i+nfit-3]/0.05d
   oplot,[-100,100],[-100,100],color=fsc_color('red')
endif else begin
   plot,tru[ok,i+ntru-3]/0.05d,fit[ok,i+nfit-3]/0.05d,xtitle=xts[i+6],ytitle=yts[i+6],$
        psym=7,charsize=2,xstyle=1,ystyle=1,$
        xrange=[-0.01d,0.01d],$
        yrange=[-0.04d,0.04d],$
        xticks=4,yticks=4
   oploterr,tru[ok,i+ntru-3]/0.05d,fit[ok,i+nfit-3]/0.05d,err[ok,i+nfit-3]/0.05d
   oplot,[-100,100],[-100,100],color=fsc_color('red')
endelse

;   if i gt 1 then begin
;      plot,tru[ok,i+ntru-3]/0.05d,(fit[ok,i+nfit-3]-tru[ok,i+ntru-3])/0.05d,$
;           xtitle=xts[i+6],ytitle=yts[i+6],$
;           psym=7,charsize=1.5,xstyle=1,ystyle=1,$
;           xrange=[-0.025d,0.025d],$
;           yrange=[-0.05d,0.05d],$
;           xticks=4,yticks=4
;      oploterr,tru[ok,i+ntru-3]/0.05d,(fit[ok,i+nfit-3]-tru[ok,i+ntru-3])/0.05d,$
;               err[ok,i+nfit-3]/0.05d
;      oplot,[-100,100],[0,0],color=fsc_color('red')
;   endif else begin
;      plot,tru[ok,i+ntru-3]/0.05d,(fit[ok,i+nfit-3]-tru[ok,i+ntru-3])/0.05d,$
;           xtitle=xts[i+6],ytitle=yts[i+6],$
;           psym=7,charsize=1.5,xstyle=1,ystyle=1,$
;           xrange=[-0.01d,0.01d],$
;           yrange=[-0.04d,0.04d],$
;           xticks=4,yticks=4
;      oploterr,tru[ok,i+ntru-3]/0.05d,(fit[ok,i+nfit-3]-tru[ok,i+ntru-3])/0.05d,$
;               err[ok,i+nfit-3]/0.05d
;      oplot,[-100,100],[0,0],color=fsc_color('red')
;   endelse

   setup_x
 ;  ps2eps,files[i+6],/delete

;   plot,fit[ok,1],(tru[ok,ntru-3+i]-fit[ok,nfit-3+i])/0.05d,xtitle=vars[6+i],/psym
;   wait,4


   use=ok[where(abs(tru[ok,ntru-3+i]-fit[ok,nfit-3+i]) lt err[ok,nfit-3+i]*2d,nuse)]

   b_m=linfit(tru[use,ntru-3+i],tru[use,ntru-3+i]-fit[use,nfit-3+i],$
              measure_errors=err[use,nfit-3+i],chisq=chisq,sigma=sig,/double)

   print,'***'
   print,vars[6+i]
   print,b_m
   print,sig
   print,chisq,'/',nok-2
   print,'***'


endfor



end
