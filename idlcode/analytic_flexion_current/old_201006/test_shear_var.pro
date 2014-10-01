pro test_shear_var, ngal, theta_e, subfrac

fieldsize=4000d

; 5 substructure radii, 10 angles
sub_radii=(fieldsize/2d)*(dindgen(5)+1d)/5d

ave_shear=dblarr(5,10)
sig_shear=dblarr(5,10)
n_used=dblarr(5,10)

; Loop over the different substructures
for k=0,4 do begin
   ; Loop over tests to get a sample
   for i=0,9 do begin
   ; Where is the substructure?
      sub_angle=2d*!dpi*randomu(seed)
      sub_ctr=fieldsize*[0.5d,0.5d]+sub_radii[k]*[cos(sub_angle),sin(sub_angle)]

   ; Galaxy positions
      gal_pos=randomu(seed,ngal,2,/double)*fieldsize

   ; Lens parameters
      main_lpars=sis_lensmodel(gal_pos,theta_e=theta_e,kappa=main_kappa,$
                               ctr=fieldsize*[0.5d,0.5d])
      sub_lpars=sis_lensmodel(gal_pos,theta_e=theta_e*subfrac,kappa=sub_kappa,$
                              ctr=sub_ctr)
      lenspars=dblarr(ngal,6)

      for j=0,ngal-1 do begin
         lenspars[j,*]=lens_superposition(main_lpars[j,*],sub_lpars[j,*],$
                                          main_kappa[j],sub_kappa[j])
      endfor

      ok=where((sqrt(total(lenspars[*,2:3]^2,2)) lt 20d) and $
               (sqrt(total(lenspars[*,4:5]^2,2)) lt 20d),nok)

      ave_shear[k,i]=mean(lenspars[ok,0:1])
      sig_shear[k,i]=stddev(lenspars[ok,0:1])
      n_used[k,i]=nok
   endfor
endfor

sig_plot_middle=0.5d*(max(sig_shear)+min(sig_shear))
sig_plot_width=(max(sig_shear)-min(sig_shear))
ave_plot_middle=0.5d*(max(ave_shear)+min(ave_shear))
ave_plot_width=(max(ave_shear)-min(ave_shear))
n_plot_middle=0.5d*(max(n_used)+min(n_used))
n_plot_width=(max(n_used)-min(n_used))

syms=[1,4,5,6,7]

plot,sub_radii[0]*(dblarr(10)+1d),ave_shear[0,*],psym=syms[0],$
     xstyle=1,xrange=[0d,0.6d]*fieldsize,ystyle=1,$
     yrange=ave_plot_middle+ave_plot_width*[-0.6d,0.6d],title='Ave g'
for i=1,4 do oplot,sub_radii[i]*(dblarr(10)+1d),ave_shear[i,*],psym=syms[i]

window,/free
plot,sub_radii[0]*(dblarr(10)+1d),sig_shear[0,*],psym=syms[0],$
     xstyle=1,xrange=[0d,0.6d]*fieldsize,ystyle=1,$
     yrange=sig_plot_middle+sig_plot_width*[-0.6d,0.6d],title='Stddev g'
for i=1,4 do oplot,sub_radii[i]*(dblarr(10)+1d),sig_shear[i,*],psym=syms[i]

window,/free
plot,sub_radii[0]*(dblarr(10)+1d),n_used[0,*],psym=syms[0],$
     xstyle=1,xrange=[0d,0.6d]*fieldsize,ystyle=1,$
     yrange=n_plot_middle+n_plot_width*[-0.6d,0.6d],title='N used'
for i=1,4 do oplot,sub_radii[i]*(dblarr(10)+1d),n_used[i,*],psym=syms[i]

end
