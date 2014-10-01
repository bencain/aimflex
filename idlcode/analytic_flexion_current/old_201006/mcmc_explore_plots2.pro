pro mcmc_explore_plots2,tag1,tag2

fit1=read_data(tag1+'_mcmc_explore_fit.dat',/quiet)
fit2=read_data(tag2+'_mcmc_explore_fit.dat',/quiet)
true1=read_data(tag1+'_mcmc_explore_true.dat',/quiet)
true2=read_data(tag2+'_mcmc_explore_true.dat',/quiet)

Gn_lim=20d

dofok1=where(fit1[*,13] gt 0d)
G1ok1=where(sqrt(total(fit1[*,8:9]^2,2)) lt Gn_lim)
G3ok1=where(sqrt(total(fit1[*,10:11]^2,2)) lt Gn_lim)
dofok2=where(fit2[*,13] gt 0d)
G1ok2=where(sqrt(total(fit2[*,8:9]^2,2)) lt Gn_lim)
G3ok2=where(sqrt(total(fit2[*,10:11]^2,2)) lt Gn_lim)

ok1=set_intersection(dofok1,G1ok1)
ok1=set_intersection(ok1,G3ok1)
ok2=set_intersection(dofok2,G1ok2)
ok2=set_intersection(ok2,G3ok2)

cs=2.5
t=10.0

plotnames=['poster'+tag1+tag2+'_g1_fixvtrue.ps',$
           'poster'+tag1+tag2+'_g2_fixvtrue.ps',$
           'poster'+tag1+tag2+'_G11_fitvtrue.ps',$
           'poster'+tag1+tag2+'_G12_fitvtrue.ps',$
           'poster'+tag1+tag2+'_G31_fitvtrue.ps',$
           'poster'+tag1+tag2+'_G32_fitvtrue.ps']

xtitles=['True g_1 values',$
         'True g_2 values',$
         'True G1_1 values (1/pixel)',$
         'True G1_2 values (1/pixel)',$
         'True G3_1 values (1/pixel)',$
         'True G3_2 values (1/pixel)']

ytitles=['Fixed g_1 values',$
         'Fixed g_2 values',$
         'Fit G1_1 values (1/pixel)',$
         'Fit G1_2 values (1/pixel)',$
         'Fit G3_1 values (1/pixel)',$
         'Fit G3_2 values (1/pixel)']

rmin=[-1d,-1d,-Gn_lim,-Gn_lim,-Gn_lim,-Gn_lim]
rmax=[1d,1d,Gn_lim,Gn_lim,Gn_lim,Gn_lim]
scale=[1d,1d,0.001d,0.001d,0.001d,0.001d]

rmin*=scale
rmax*=scale

for i=0,5 do begin
   setup_ps,plotnames[i]
   plot,true1[ok1,i+6]*scale[i],fit1[ok1,i+6]*scale[i],$
        xtitle=xtitles[i],ytitle=ytitles[i],$
        psym=5,xrange=[rmin[i],rmax[i]],yrange=[rmin[i],rmax[i]],$
        xstyle=1,ystyle=1
   oplot,true2[ok2,i+6]*scale[i],fit2[ok2,i+6]*scale[i],$
         psym=7
   oplot,[rmin[i],rmax[i]],[rmin[i],rmax[i]]
   xyouts,0.2d*rmax[i],0.8*rmin[i],'sigma_tri = '+$
          strcompress(string(stddev(true1[ok1,i+6]*scale[i]-$
                                    fit1[ok1,i+6]*scale[i]),$
                             format='(e10.2)'),/remove_all)
   xyouts,0.2d*rmax[i],0.9*rmin[i],'sigma_x = '+$
          strcompress(string(stddev(true2[ok2,i+6]*scale[i]-$
                                    fit2[ok2,i+6]*scale[i]),$
                             format='(e10.2)'),/remove_all)
   
   setup_x
   ps2eps,plotnames[i],/delete
endfor

end
