pro phys7a_plot

eps=3d
sig=4d
r=3*sig*(1d + dindgen(1000))/1000d
v=4d*eps*((sig/r)^12 - (sig/r)^6)

set_plot,'PS'
device,filename='Q1fig.ps';,xsize=7,ysize=5,xoffset=0.5,yoffset=1.0,/inches

;plot,r,v,xrange=[0,12d],yrange=[-5d,5d],xtitle='Separation (10!S!E-10!R  m)',ytitle='Potential Energy (10!S!E-21!R  J)',charsize=2,charthick=3,thick=5,xstyle=1,ystyle=1,xgridstyle=1,ygridstyle=1,ticklen=1d

e=dindgen(21)
t=dindgen(16)

plot,e,t,xrange=[0,20d],yrange=[0d,15d],xstyle=1,ystyle=1,ticklen=1d,/isotropic,psym=3,$
     xticks=20,yticks=15,xtickname=REPLICATE(' ', 21),ytickname=REPLICATE(' ', 16)

setup_x


end
