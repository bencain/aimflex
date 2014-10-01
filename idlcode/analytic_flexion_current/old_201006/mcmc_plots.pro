pro MCMC_PLOTS,tag,wsfitpars,nsfitpars,truepars,pnames


; Now we do the plots.  We want to make the following 8 plots:
; E1-alpha
; E1-g1
; E1-g2
; E2-g2
; E2-g1
; E2-alpha
; g1-alpha
; g2-alpha

nplot=8
firstpars =[4,4,4,5,5,5,6,7]
secondpars=[3,6,7,7,6,3,3,3]

; Plot the fits
for i=0,nplot-1 do begin

   forloop_status,i,nplot,wnum,label='Plotting par vs par to .eps files...'
;   print,'plot',i+1

   wsfile=tag+'wsplot'+strcompress(string(i+1),/remove_all)+'.ps'
   nsfile=tag+'nsplot'+strcompress(string(i+1),/remove_all)+'.ps'

   setup_ps,wsfile

   plot,wsfitpars[*,firstpars[i]],wsfitpars[*,secondpars[i]],$
        xtitle=pnames[firstpars[i]],ytitle=pnames[secondpars[i]],$
        psym=7, charsize=3.0,xmargin=[30,15], ymargin=[15,8],$
        title=$
        pnames[firstpars[i]]+' vs '+pnames[secondpars[i]]+' - with shear',$
        charthick=3.0, thick=3.0, xthick=3.0, ythick=3.0

   setup_x
   ps2eps,wsfile,/delete

   setup_ps,nsfile
   
   plot,nsfitpars[*,firstpars[i]],nsfitpars[*,secondpars[i]],$
        xtitle=pnames[firstpars[i]],ytitle=pnames[secondpars[i]],$
        psym=7, charsize=3.0,xmargin=[30,15], ymargin=[15,8],$
        title=$
        pnames[firstpars[i]]+' vs '+pnames[secondpars[i]]+' - no shear',$
        charthick=3.0, thick=3.0, xthick=3.0, ythick=3.0

   setup_x
   ps2eps,nsfile,/delete

endfor

forloop_status,0,0,wnum,/delete

; Plot fits vs true values for each parameter
for i=0,11 do begin
   
   forloop_status,i,12,wnum,label='Plotting fit vs true to .eps files...'

;   print,'fvt',i+1

   if i+1 lt 10 then gap='0' else gap=''

   wsfile=tag+'wsfvtplot'+gap+strcompress(string(i+1),/remove_all)+'.ps'
   nsfile=tag+'nsfvtplot'+gap+strcompress(string(i+1),/remove_all)+'.ps'

   setup_ps,wsfile

   xmin=min(truepars[*,i])
   xmax=max(truepars[*,i])
   range=xmax-xmin
   midpt=(xmax+xmin)/2d
   pts=midpt+1.5d*[-range/2d,range/2d]

   if pts[0] eq pts[1] then begin
      if pts[0] ne 0d then pts=pts[0]*[0.5d,1.5d] else pts=[-1d-6,1d-6]
   endif

   style=1

;   print,pts

   plot,truepars[*,i],wsfitpars[*,i],$
        xtitle='True '+pnames[i],ytitle='Fit '+pnames[i],$
        psym=7, charsize=3.0,xmargin=[30,30], ymargin=[15,15],$
        title='Fit vs True - With shear - '+pnames[i],$
        charthick=3.0, thick=3.0, xthick=3.0, ythick=3.0,$
        xrange=pts,xstyle=style
   oplot,pts,pts

   setup_x
   ps2eps,wsfile,/delete

   setup_ps,nsfile

   plot,truepars[*,i],nsfitpars[*,i],$
        xtitle='True '+pnames[i],ytitle='Fit '+pnames[i],$
        psym=7, charsize=3.0,xmargin=[30,30], ymargin=[15,15],$
        title='Fit vs True - No shear - '+pnames[i],$
        charthick=3.0, thick=3.0, xthick=3.0, ythick=3.0,$
        xrange=pts,xstyle=style
   oplot,pts,pts

   setup_x
   ps2eps,nsfile,/delete

endfor

forloop_status,0,0,wnum,/delete



; Plot chi-squared histograms

setup_ps,tag+'nshistplot.ps'

plot_hist,nsfitpars[*,12],xtitle='Chi^2',ytitle='Fraction',$
          charsize=3.0,xmargin=[30,30], ymargin=[15,15],$
          charthick=3.0, thick=3.0, xthick=3.0, ythick=3.0,$
          title='Chi-squared - No Lensing',$
          ystyle=1,/norm

setup_x
ps2eps,tag+'nshistplot.ps',/delete

setup_ps,tag+'wshistplot.ps'

plot_hist,wsfitpars[*,12],xtitle='Chi^2',ytitle='Fraction',$
          charsize=3.0,xmargin=[30,30], ymargin=[15,15],$
          charthick=3.0, thick=3.0, xthick=3.0, ythick=3.0,$
          title='Chi-squared - With Lensing',$
          ystyle=1,/norm

setup_x
ps2eps,tag+'wshistplot.ps',/delete


; Plot chi-squared vs iteration number

setup_ps,tag+'ns_chisqvniter_plot.ps'

plot,nsfitpars[*,13],nsfitpars[*,12],ytitle='Chi-Squared',xtitle='Niter',$
     charsize=3.0,xmargin=[30,30], ymargin=[15,15],$
     charthick=3.0, thick=3.0, xthick=3.0, ythick=3.0,$
     title='Chi^2 vs N_iteration - No Lensing',psym=7,/ylog

setup_x
ps2eps,tag+'ns_chisqvniter_plot.ps',/delete

setup_ps,tag+'ws_chisqvniter_plot.ps'

plot,wsfitpars[*,13],wsfitpars[*,12],ytitle='Chi-Squared',xtitle='Niter',$
     charsize=3.0,xmargin=[30,30], ymargin=[15,15],$
     charthick=3.0, thick=3.0, xthick=3.0, ythick=3.0,$
     title='Chi^2 vs N_iteration - With Lensing',psym=7,/ylog

setup_x
ps2eps,tag+'ws_chisqvniter_plot.ps',/delete

end
