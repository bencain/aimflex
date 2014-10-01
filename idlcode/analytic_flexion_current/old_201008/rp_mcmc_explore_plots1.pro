pro rp_mcmc_explore_plots1, tag

; Filenames from RP_MCMC_EXPLORE.PRO
startfile=tag+'_start.dat'
truefile= tag+'_true.dat'
errfile = tag+'_errs.dat'
fitfile = tag+'_fitps.dat'
fomfile = tag+'_foms.dat'

startps=read_data(startfile,/quiet)
trueps=read_data(truefile,/quiet)
errs=read_data(errfile,/quiet)
foms=read_data(fomfile,/quiet)
fitps=read_data(fitfile,/quiet,size=fitsize)

; We'll cut out the fits that failed, ran to the end of the
; allowed iterations, or pegged a parameter.
notfailed=where(foms[*,8] eq 0, n_notfailed)
notitered=where(foms[*,7] lt 250,n_notitered)
notpegged=where(foms[*,6] eq 0,n_notpegged)

ok=set_intersection(notfailed,notitered)
ok=set_intersection(ok,notpegged,count=ngood)

print,'  '
print,'Failed fits:.................'+string(fitsize[0]-n_notfailed)
print,'Fits with pegged parameters:.'+string(fitsize[0]-n_notpegged)
print,'Maxed out on iterations :....'+string(fitsize[0]-n_notitered)
print,'Good fits:...................'+string(ngood)
print,'  '

; Make plots
;  chi^2/DoF vs each of the parameters (12)
;  Fit vs true values for each parameter (12)
;
psfiles=['_chisq_vs_logN0.ps',$
          '_chisq_vs_Xc.ps',$
          '_chisq_vs_Yc.ps',$
          '_chisq_vs_alpha.ps',$
          '_chisq_vs_ecross.ps',$
          '_chisq_vs_eplus.ps',$ 
          '_chisq_vs_g1.ps',$
          '_chisq_vs_g2.ps',$
          '_chisq_vs_G11.ps',$
          '_chisq_vs_G12.ps',$
          '_chisq_vs_G31.ps',$
          '_chisq_vs_G32.ps',$
          '_fitlogN0_vs_truelogN0.ps',$
          '_fitXc_vs_trueXc.ps',$
          '_fitYc_vs_trueYc.ps',$
          '_fitalpha_vs_truealpha.ps',$
          '_fiteplus_vs_trueeplus.ps',$
          '_fitecross_vs_trueecross.ps',$
          '_fitg1_vs_trueg1.ps',$
          '_fitg2_vs_trueg2.ps',$
          '_fitG11_vs_trueG11.ps',$
          '_fitG12_vs_trueG12.ps',$
          '_fitG31_vs_trueG31.ps',$
          '_fitG32_vs_trueG32.ps']
psfiles=tag+psfiles

xtitles=['True logN0 (counts)',$
         'True Xc (pixels)',$
         'True Yc (pixels)',$
         'True alpha (pixels)',$
         'True E+',$
         'True Ex',$         
         'True g_1',$
         'True g_2',$
         'True G_11 (1/pixel)',$
         'True G_12 (1/pixel)',$
         'True G_31 (1/pixel)',$
         'True G_32 (1/pixel)',$
         'True logN0 (counts)',$
         'True Xc (pixels)',$
         'True Yc (pixels)',$
         'True alpha (pixels)',$
         'True E+',$
         'True Ex',$
         'True g_1',$
         'True g_2',$
         'True G_11 (1/pixel)',$
         'True G_12 (1/pixel)',$
         'True G_31 (1/pixel)',$
         'True G_32 (1/pixel)']

ytitles=['Chi-squared per Degree of Freedom',$
         'Chi-squared per Degree of Freedom',$
         'Chi-squared per Degree of Freedom',$
         'Chi-squared per Degree of Freedom',$
         'Chi-squared per Degree of Freedom',$
         'Chi-squared per Degree of Freedom',$
         'Chi-squared per Degree of Freedom',$
         'Chi-squared per Degree of Freedom',$
         'Chi-squared per Degree of Freedom',$
         'Chi-squared per Degree of Freedom',$
         'Chi-squared per Degree of Freedom',$
         'Chi-squared per Degree of Freedom',$
         'Fit logN0 (counts)',$
         'Fit Xc (pixels)',$
         'Fit Yc (pixels)',$
         'Fit alpha (pixels)',$
         'Fit E+',$
         'Fit Ex',$
         'Fit g_1',$
         'Fit g_2',$
         'Fit G_11 (1/pixel)',$
         'Fit G_12 (1/pixel)',$
         'Fit G_31 (1/pixel)',$
         'Fit G_32 (1/pixel)']

xs=dblarr(n_elements(xtitles),ngood)
ys=xs

; True pars vs chi^2/Dof, Fit vs True pars for the first 24 plots. Put
; in the xs and ys for these
for i=0,11 do begin
   xs[i,*]=trueps[ok,i+1]
   xs[i+12,*]=trueps[ok,i+1]
endfor

for i=0,11 do begin
   ys[i,*]=foms[ok,1]/foms[ok,2]
   ys[i+12,*]=fitps[ok,i+1]
endfor

for i=0,n_elements(psfiles)-1 do begin
   setup_ps,psfiles[i]
   plot,xs[i,*],ys[i,*],psym=1,ylog=(i lt 12),$
        xtitle=xtitles[i],ytitle=ytitles[i]
   if i gt 11 then oplot,[min(xs[i,*]),max(xs[i,*])],$
                        [min(xs[i,*]),max(xs[i,*])]
   setup_x
endfor

for i=0,n_elements(psfiles)-1 do ps2eps,psfiles[i],/delete

end
