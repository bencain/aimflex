pro todo2

dir='data/paper2011/'

fit=read_data(dir+'run02_fit.dat',/quiet)
tru=read_data(dir+'run02_true.dat',/quiet)
fom=read_data(dir+'run02_fom.dat',/quiet)
err=read_data(dir+'run02_err.dat',/quiet)

ok=where((fom[*,3] lt 1.5d) and $ ; Chisq/dof ok
         (fom[*,11] eq 0) and $   ; no flags
         (fit[*,1] gt 0d) and $   ; flux limit
         (fit[*,4] gt 1.8d) and $ ; size limit
         (err[*, 9] gt 1d-4) and $ ; error limits
         (err[*,10] gt 1d-4) and $
         (err[*,11] gt 1d-4) and $
         (err[*,12] gt 1d-4) and $
         (err[*, 9] lt 1d-3) and $
         (err[*,10] lt 1d-3) and $
         (err[*,11] lt 1d-3) and $
         (err[*,12] lt 1d-3),nok)

print, nok

sis=lenspars2sispars(fit[ok,7:12])

; delta(psi_mn) vs sig(psi_mn)
plot, err[*,9:12],abs(fit[*,9:12]-tru[*,9:12]),/xlog,/ylog,/psym
oplot,[1d-10,1d5],[1d-10,1d5]

; delta(psi_1n) vs abs(psi_1n)
;plot, tru[*,4], sqrt(total((fit[*,9:10]-tru[*,9:10])^2,2)),/psym,/ylog


;plot, sis[*,0],abs(fit[ok, 9]-tru[ok, 9])/fom[ok,6],/psym,/ylog
;oplot,sis[*,0],abs(fit[ok,10]-tru[ok,10])/fom[ok,6],/psym
;oplot,sis[*,0],abs(fit[ok,11]-tru[ok,11])/fom[ok,6],/psym
;oplot,sis[*,0],abs(fit[ok,12]-tru[ok,12])/fom[ok,6],/psym

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

dir='data/paper2011/'
n=1000L
sig_I=1.1d-3
sig_g=0.2d

;aim_mcmc_explore, n, sig_I, dir+'run01', /const_error, cutrad_scale=1.1d, shearnoise=sig_g, /sersic_true
;aim_mcmc_explore, n, sig_I, dir+'run02', /const_error, cutrad_scale=1.1d, shearnoise=sig_g, /sersic_true, /sis_lens

head1=dir+'run01_'
head2=dir+'run02_'
aim_mcmc_plots,head1+'fit.dat',$
               head1+'fom.dat',$
               head1+'true.dat',$
               head1+'err.dat',$
               head1+'start.dat',$
               tag=head1+'plot'
aim_mcmc_plots,head2+'fit.dat',$
               head2+'fom.dat',$
               head2+'true.dat',$
               head2+'err.dat',$
               head2+'start.dat',$
               tag=head2+'plot'

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


dir='data/a1689_hla/fourfilter/aim_fit/'


fit=read_data(dir+'use_fit.dat',/quiet)
fom=read_data(dir+'use_fom.dat',/quiet)
sta=read_data(dir+'use_start.dat',/quiet)
err=read_data(dir+'use_err.dat',/quiet)
cor=read_data(dir+'use_corr.dat',/quiet)

ok=where((err[*, 9] gt 1d-4) and (err[*,10] gt 1d-4) and $
         (err[*,11] gt 1d-4) and (err[*,12] gt 1d-4),nok)

fit=fit[ok,*]
fom=fom[ok,*]
sta=sta[ok,*]
err=err[ok,*]
cor=cor[ok,*]


;tag='img_'
;E=sqrt(total(sta[*,5:6]^2,2))
;eps=(1d - E)/(1d + E)
;a=sta[*,4]/sqrt(eps)

tag='src_'
E=sqrt(total(fit[*,7:8]^2,2))
eps=(1d - E)/(1d + E)
a=fit[*,6]/sqrt(eps)


fflex=dcomplex(fit[*,11],fit[*,12])
gflex=dcomplex(fit[*,13],fit[*,14])

ctr=dcomplex(3040d,2655d)
pos=dcomplex(fit[*,1],fit[*,2]) - ctr
rad=abs(pos)
ph=pos/rad

dF=dblarr(nok)
dG=dblarr(nok)
for i=0,nok-1 do begin
   dF[i]=max(err[i,9:10])
   dG[i]=max(err[i,11:12])
endfor

lo_Ferr=(sort(dF))[0:floor(nok/2d)-1]
hi_Ferr=(sort(dF))[floor(nok/2d):nok-1]
print,'F Err',max(dF[lo_Ferr])/0.05d

lo_Gerr=(sort(dG))[0:floor(nok/2d)-1]
hi_Gerr=(sort(dG))[floor(nok/2d):nok-1]
print,'G Err',max(dG[lo_Gerr])/0.05d

lo_rad=(sort(rad))[0:floor(nok/2d)-1]
hi_rad=(sort(rad))[floor(nok/2d):nok-1]

print,'Radius',max(rad[lo_rad])*0.05d

;;;;;;;;;;;; 1 - FLEXION
;;;;;;;;;;;;
;;;;;;;;;;;;; Error splitting ;;;;;;;;;;;;;;;;;;;;;;;

x=[a*real_part(fflex),a*imaginary(fflex)]
ave=mean(x)
sig=stddev(x)
setup_ps,tag+'aFn_all.ps'
plot_hist,x,xtitle='a*F_n',xrange=[-0.2,0.2],yrange=[0,50],binsize=0.005,charsize=1.5
xyouts,0.75,0.75,'mean = '+strcompress(string(ave,format='(F7.3)'),/remove_all),/norm
xyouts,0.75,0.70,'stddev = '+strcompress(string(sig,format='(F7.3)'),/remove_all),/norm
setup_x
ps2eps,tag+'aFn_all.ps',/delete


x=[(a*real_part(fflex))[lo_Ferr],(a*imaginary(fflex))[lo_Ferr]]
xave=mean(x)
xsig=stddev(x)
y=[(a*real_part(fflex))[hi_Ferr],(a*imaginary(fflex))[hi_Ferr]]
yave=mean(y)
ysig=stddev(y)
setup_ps,tag+'aFn_err.ps'
plot_hist,x,xtitle='a*F_n',xrange=[-0.2,0.2],yrange=[0,30],binsize=0.005,charsize=1.5
plot_hist,y,binsize=0.005,/overplot,color=fsc_color('blue')
xyouts,0.75,0.85,'Low sig(F_n)',/norm
xyouts,0.75,0.80,'mean = '+strcompress(string(xave,format='(F7.3)'),/remove_all),/norm
xyouts,0.75,0.75,'stddev = '+strcompress(string(xsig,format='(F7.3)'),/remove_all),/norm

xyouts,0.75,0.70,'High sig(F_n)',/norm,color=fsc_color('blue')
xyouts,0.75,0.65,'mean = '+strcompress(string(yave,format='(F7.3)'),/remove_all),$
       /norm,color=fsc_color('blue')
xyouts,0.75,0.60,'stddev = '+strcompress(string(ysig,format='(F7.3)'),/remove_all),$
       /norm,color=fsc_color('blue')
setup_x
ps2eps,tag+'aFn_err.ps',/delete


;;;;;;;;;;;;; Radial stuff ::::::::::::::::::::

x=a*real_part(fflex/ph)
ave=mean(x)
sig=stddev(x)
setup_ps,tag+'aFr_all.ps'
plot_hist,x,xtitle='a*F_r',xrange=[-0.2,0.2],yrange=[0,30],binsize=0.005,charsize=1.5
xyouts,0.75,0.75,'mean = '+strcompress(string(ave,format='(F7.3)'),/remove_all),/norm
xyouts,0.75,0.70,'stddev = '+strcompress(string(sig,format='(F7.3)'),/remove_all),/norm
setup_x
ps2eps,tag+'aFr_all.ps',/delete

x=(a*real_part(fflex/ph))[lo_rad]
xave=mean(x)
xsig=stddev(x)
y=(a*real_part(fflex/ph))[hi_rad]
yave=mean(y)
ysig=stddev(y)
setup_ps,tag+'aFr_rad.ps'
plot_hist,x,xtitle='a*F_r',xrange=[-0.2,0.2],yrange=[0,16],binsize=0.005,charsize=1.5
plot_hist,y,binsize=0.005,/overplot,color=fsc_color('blue')
xyouts,0.75,0.85,'Low R',/norm
xyouts,0.75,0.80,'mean = '+strcompress(string(xave,format='(F7.3)'),/remove_all),/norm
xyouts,0.75,0.75,'stddev = '+strcompress(string(xsig,format='(F7.3)'),/remove_all),/norm

xyouts,0.75,0.70,'High R',/norm,color=fsc_color('blue')
xyouts,0.75,0.65,'mean = '+strcompress(string(yave,format='(F7.3)'),/remove_all),$
       /norm,color=fsc_color('blue')
xyouts,0.75,0.60,'stddev = '+strcompress(string(ysig,format='(F7.3)'),/remove_all),$
       /norm,color=fsc_color('blue')
setup_x
ps2eps,tag+'aFr_rad.ps',/delete


x=a*imaginary(fflex/ph)
ave=mean(x)
sig=stddev(x)
setup_ps,tag+'aFt_all.ps'
plot_hist,x,xtitle='a*F_t',xrange=[-0.2,0.2],yrange=[0,30],binsize=0.005,charsize=1.5
xyouts,0.75,0.75,'mean = '+strcompress(string(ave,format='(F7.3)'),/remove_all),/norm
xyouts,0.75,0.70,'stddev = '+strcompress(string(sig,format='(F7.3)'),/remove_all),/norm
setup_x
ps2eps,tag+'aFt_all.ps',/delete

x=(a*imaginary(fflex/ph))[lo_rad]
xave=mean(x)
xsig=stddev(x)
y=(a*imaginary(fflex/ph))[hi_rad]
yave=mean(y)
ysig=stddev(y)
setup_ps,tag+'aFt_rad.ps'
plot_hist,x,xtitle='a*F_t',xrange=[-0.2,0.2],yrange=[0,16],binsize=0.005,charsize=1.5
plot_hist,y,binsize=0.005,/overplot,color=fsc_color('blue')
xyouts,0.75,0.85,'Low R',/norm
xyouts,0.75,0.80,'mean = '+strcompress(string(xave,format='(F7.3)'),/remove_all),/norm
xyouts,0.75,0.75,'stddev = '+strcompress(string(xsig,format='(F7.3)'),/remove_all),/norm

xyouts,0.75,0.70,'High R',/norm,color=fsc_color('blue')
xyouts,0.75,0.65,'mean = '+strcompress(string(yave,format='(F7.3)'),/remove_all),$
       /norm,color=fsc_color('blue')
xyouts,0.75,0.60,'stddev = '+strcompress(string(ysig,format='(F7.3)'),/remove_all),$
       /norm,color=fsc_color('blue')
setup_x
ps2eps,tag+'aFt_rad.ps',/delete



;;;;; Magnitude stuff

x=a*abs(fflex)
setup_ps,tag+'aF_abs.ps'
plot_hist,x,xtitle='a*|F|',xrange=[0,0.2],yrange=[0,30],binsize=0.005,charsize=1.5
setup_x
ps2eps,tag+'aF_abs.ps',/delete

;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;; G-flexion
;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;; 1 - FLEXION
;;;;;;;;;;;;
;;;;;;;;;;;;; Error splitting ;;;;;;;;;;;;;;;;;;;;;;;

x=[a*real_part(gflex),a*imaginary(gflex)]
ave=mean(x)
sig=stddev(x)
setup_ps,tag+'aGn_all.ps'
plot_hist,x,xtitle='a*G_n',xrange=[-0.2,0.2],yrange=[0,50],binsize=0.005,charsize=1.5
xyouts,0.75,0.75,'mean = '+strcompress(string(ave,format='(F7.3)'),/remove_all),/norm
xyouts,0.75,0.70,'stddev = '+strcompress(string(sig,format='(F7.3)'),/remove_all),/norm
setup_x
ps2eps,tag+'aGn_all.ps',/delete


x=[(a*real_part(gflex))[lo_Gerr],(a*imaginary(gflex))[lo_Gerr]]
xave=mean(x)
xsig=stddev(x)
y=[(a*real_part(gflex))[hi_Gerr],(a*imaginary(gflex))[hi_Gerr]]
yave=mean(y)
ysig=stddev(y)
setup_ps,tag+'aGn_err.ps'
plot_hist,x,xtitle='a*G_n',xrange=[-0.2,0.2],yrange=[0,30],binsize=0.005,charsize=1.5
plot_hist,y,binsize=0.005,/overplot,color=fsc_color('blue')
xyouts,0.75,0.85,'Low sig(G_n)',/norm
xyouts,0.75,0.80,'mean = '+strcompress(string(xave,format='(F7.3)'),/remove_all),/norm
xyouts,0.75,0.75,'stddev = '+strcompress(string(xsig,format='(F7.3)'),/remove_all),/norm

xyouts,0.75,0.70,'High sig(G_n)',/norm,color=fsc_color('blue')
xyouts,0.75,0.65,'mean = '+strcompress(string(yave,format='(F7.3)'),/remove_all),$
       /norm,color=fsc_color('blue')
xyouts,0.75,0.60,'stddev = '+strcompress(string(ysig,format='(F7.3)'),/remove_all),$
       /norm,color=fsc_color('blue')
setup_x
ps2eps,tag+'aGn_err.ps',/delete


;;;;;;;;;;;;; Radial stuff ::::::::::::::::::::

x=a*real_part(gflex/ph)
ave=mean(x)
sig=stddev(x)
setup_ps,tag+'aGr_all.ps'
plot_hist,x,xtitle='a*G_r',xrange=[-0.2,0.2],yrange=[0,30],binsize=0.005,charsize=1.5
xyouts,0.75,0.75,'mean = '+strcompress(string(ave,format='(F7.3)'),/remove_all),/norm
xyouts,0.75,0.70,'stddev = '+strcompress(string(sig,format='(F7.3)'),/remove_all),/norm
setup_x
ps2eps,tag+'aGr_all.ps',/delete

x=(a*real_part(gflex/ph))[lo_rad]
xave=mean(x)
xsig=stddev(x)
y=(a*real_part(gflex/ph))[hi_rad]
yave=mean(y)
ysig=stddev(y)
setup_ps,tag+'aGr_rad.ps'
plot_hist,x,xtitle='a*G_r',xrange=[-0.2,0.2],yrange=[0,16],binsize=0.005,charsize=1.5
plot_hist,y,binsize=0.005,/overplot,color=fsc_color('blue')
xyouts,0.75,0.85,'Low R',/norm
xyouts,0.75,0.80,'mean = '+strcompress(string(xave,format='(F7.3)'),/remove_all),/norm
xyouts,0.75,0.75,'stddev = '+strcompress(string(xsig,format='(F7.3)'),/remove_all),/norm

xyouts,0.75,0.70,'High R',/norm,color=fsc_color('blue')
xyouts,0.75,0.65,'mean = '+strcompress(string(yave,format='(F7.3)'),/remove_all),$
       /norm,color=fsc_color('blue')
xyouts,0.75,0.60,'stddev = '+strcompress(string(ysig,format='(F7.3)'),/remove_all),$
       /norm,color=fsc_color('blue')
setup_x
ps2eps,tag+'aGr_rad.ps',/delete


x=a*imaginary(gflex/ph)
ave=mean(x)
sig=stddev(x)
setup_ps,tag+'aGt_all.ps'
plot_hist,x,xtitle='a*G_t',xrange=[-0.2,0.2],yrange=[0,30],binsize=0.005,charsize=1.5
xyouts,0.75,0.75,'mean = '+strcompress(string(ave,format='(F7.3)'),/remove_all),/norm
xyouts,0.75,0.70,'stddev = '+strcompress(string(sig,format='(F7.3)'),/remove_all),/norm
setup_x
ps2eps,tag+'aGt_all.ps',/delete

x=(a*imaginary(gflex/ph))[lo_rad]
xave=mean(x)
xsig=stddev(x)
y=(a*imaginary(gflex/ph))[hi_rad]
yave=mean(y)
ysig=stddev(y)
setup_ps,tag+'aGt_rad.ps'
plot_hist,x,xtitle='a*G_t',xrange=[-0.2,0.2],yrange=[0,16],binsize=0.005,charsize=1.5
plot_hist,y,binsize=0.005,/overplot,color=fsc_color('blue')
xyouts,0.75,0.85,'Low R',/norm
xyouts,0.75,0.80,'mean = '+strcompress(string(xave,format='(F7.3)'),/remove_all),/norm
xyouts,0.75,0.75,'stddev = '+strcompress(string(xsig,format='(F7.3)'),/remove_all),/norm

xyouts,0.75,0.70,'High R',/norm,color=fsc_color('blue')
xyouts,0.75,0.65,'mean = '+strcompress(string(yave,format='(F7.3)'),/remove_all),$
       /norm,color=fsc_color('blue')
xyouts,0.75,0.60,'stddev = '+strcompress(string(ysig,format='(F7.3)'),/remove_all),$
       /norm,color=fsc_color('blue')
setup_x
ps2eps,tag+'aGt_rad.ps',/delete



;;;;; Magnitude stuff

x=a*abs(gflex)
setup_ps,tag+'aG_abs.ps'
plot_hist,x,xtitle='a*|G|',xrange=[0,0.2],yrange=[0,30],binsize=0.005,charsize=1.5
setup_x
ps2eps,tag+'aG_abs.ps',/delete



end
