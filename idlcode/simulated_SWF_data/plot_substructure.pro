pro plot_substructure

te=30d ; Main halo Einstein radius in arcsec

; In units of the main halo Einstein radius
r = te*2d*(dindgen(1d3)+1.0000120123d)/1d3


; Substructure at 1.4 main-halo Einstein radii
rsub=te*0.90000912312932131d
; Einstein radius of the substructure as a fraction of the main halo
tsub=0.75d*te/30d

f_halo=-te*0.5d/r^2
k_halo=te*0.5d/r

f_sub=-0.5d*tsub/((r-rsub)*abs(r-rsub))
k_sub=0.5*tsub/abs(r-rsub)

f_tot=f_halo+f_sub
k_tot=k_halo+k_sub

psi_halo=f_halo/(1d - k_halo)
psi_sub=f_sub/(1d - k_sub)
psi_tot=f_tot/(1d - k_tot)

below=where(r lt rsub)
above=where(r gt rsub)
nsub=5d

a=0.5d

set_plot,'PS'
!p.font=0
device,filename='substr.eps',xsize=6,ysize=6,/inches,decomposed=1, color=1, bits_per_pixel=8,/encapsulated
plot,r[below],a*(f_tot[below]),xrange=[5,40],yrange=[-10,10]/te, thick=5,xstyle=1,ystyle=1,$
     xtitle='r (arcsec)',ytitle='Scaled Flexion a!8F!3',charsize=1.5,charthick=4
oplot,r[above],a*(f_tot[above]), thick=5
oplot,r,a*(f_halo),linestyle=2, thick=5, color=fsc_color('red')
oploterror,19.5d,0.15d,0.03d,0.03d, errthick=5
xyouts,8d,0.14d,'Intrinsic Scatter',charthick=4,charsize=1.1
device,/close
set_plot,'X'



b=30d
s=1d
q=0.5d
a=!dpi/2d

theta=dblarr(1d4,2)
theta[*,0]=0.1d*b*dindgen(1d4)/1d4

lens=[0d,0d,b,s,q,a]
fields=nie_lens(theta,lens)
plot,theta[*,0],sqrt(total(fields[*,4:5]^2,2))

for i=0,4 do begin
   a-=!dpi/10d
   print,a
   lens=[0d,0d,b,s,q,a]
   fields=nie_lens(theta,lens)
   oplot,theta[*,0],sqrt(total(fields[*,4:5]^2,2))
endfor

end
