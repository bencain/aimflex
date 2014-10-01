pro mcmc_explore_plots1, tag, Gn_lim=Gn_lim, flexionscale=flexionscale

startfile=tag+'_mcmc_explore_start.dat'
truefile=tag+'_mcmc_explore_true.dat'
errfile=tag+'_mcmc_explore_err.dat'
fitfile=tag+'_mcmc_explore_fit.dat'

startps=read_data(startfile,/quiet)
trueps=read_data(truefile,/quiet)
errs=read_data(errfile,/quiet)
fitps=read_data(fitfile,/quiet,size=fitsize)


; Make several plots
;chisq/dof vs I0, alpha, and lensing parameters

if not keyword_set(Gn_lim) then Gn_lim=50d
if not keyword_set(flexionscale) then flexionscale=1d3

dofok=where(fitps[*,13] gt 0d,ncomplement=nbaddof)
G1ok=where(sqrt(total(fitps[*,8:9]^2,2)) lt Gn_lim,ncomplement=nbadG1)
G3ok=where(sqrt(total(fitps[*,10:11]^2,2)) lt Gn_lim,ncomplement=nbadG3)

ok=set_intersection(dofok,G1ok)
ok=set_intersection(ok,G3ok,count=ngood)

print,strcompress(string(nbaddof),/remove_all)+' have negative DoF.'
print,strcompress(string(nbadG1),/remove_all)+' have |G1| > '+$
      strcompress(string(Gn_lim),/remove_all)
print,strcompress(string(nbadG3),/remove_all)+' have |G3| > '+$
      strcompress(string(Gn_lim),/remove_all)
print,strcompress(string(ngood),/remove_all)+' are okay.'

psfiles=strarr(24)+tag
psfiles+=['_chisq_vs_I0.ps',$
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
          '_fitI0_vs_trueI0.ps',$
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

xtitles=['True I0 (counts/square pixel)',$
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
         'True I0 (counts/square pixel)',$
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
         'Fit I0 (counts/square pixel)',$
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

for i=0,11 do begin
   if i gt 7 then scale=flexionscale else scale=1d
   xs[i,*]=trueps[ok,i]/scale
   xs[i+12,*]=trueps[ok,i]/scale
endfor

for i=0,11 do ys[i,*]=fitps[ok,12]/fitps[ok,13]
for i=0,11 do begin
   if i gt 7 then scale=flexionscale else scale=1d
   ys[i+12,*]=fitps[ok,i]/scale
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
