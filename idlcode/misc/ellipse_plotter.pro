pro ellipse_plotter, n_ellipse, tag

;!p.multi=[0,1,2]
;window,xsize=450,ysize=900

; Ellipse Center Positions
size=1000d
xpos=randomu(seed,n_ellipse)*size
ypos=randomu(seed,n_ellipse)*size

;xpos=(0.5d +(dindgen(n_ellipse) mod 5))*size/5d
;ypos=(0.5d +double(lindgen(n_ellipse)/5L))*size/5d

; Ellipse Size, Axis Ratio and Position Angles
r0=20d
ecc=0.7d
pas=!dpi*randomu(seed,n_ellipse)-!dpi/2d

phis=2d*!dpi*dindgen(300)/299d

; Lens stuff
theta_e=size/5d
r_lens=sqrt((xpos-size/2d)^2+(ypos-size/2d)^2)
ang_lens=atan(ypos-size/2d,xpos-size/2d)

gs=-theta_e*dcomplex(cos(2d*ang_lens),sin(2d*ang_lens))/(2d*r_lens-theta_e)

chi_arr=dblarr(n_ellipse,4)
chi_arr[*,0]=xpos
chi_arr[*,1]=ypos

gs_arr=chi_arr
gs_arr[*,2]=real_part(gs)
gs_arr[*,3]=imaginary(gs)

; Plot w/o lensing
setup_ps,tag+'_orig.ps'
plot,[0,size,size,0,0],[0,0,size,size,0],xstyle=4,ystyle=4,/isotropic
for i=0,n_ellipse-1 do begin
   oplot,$
      xpos[i]+r0*(1d -ecc^2)*cos(phis)/(1d +ecc*cos(phis-pas[i]))+$
      r0*ecc*cos(pas[i]),$
      ypos[i]+r0*(1d -ecc^2)*sin(phis)/(1d +ecc*cos(phis-pas[i]))+$
      r0*ecc*sin(pas[i]), thick=5d
endfor
setup_x

; Plot w/ lensing
setup_ps,tag+'_lens.ps'
plot,[0,size,size,0,0],[0,0,size,size,0],xstyle=4,ystyle=4,/isotropic
for i=0,n_ellipse-1 do begin
   
   eps=sqrt(1d -ecc^2)
   chi_src=(1d -eps^2)*dcomplex(cos(2d*pas[i]),2d*sin(pas[i]))/(1d +eps^2)

   chi=(chi_src+2d*gs[i]+gs[i]^2*conj(chi_src))/$
       (1d +gs[i]*conj(gs[i])+2d*real_part(gs[i]*conj(chi_src)))

   chi_arr[i,2]=real_part(chi)
   chi_arr[i,3]=imaginary(chi)

   lpa=0.5d*atan(imaginary(chi),real_part(chi))
   leps=sqrt((1d - sqrt(imaginary(chi)^2+real_part(chi)^2))/$
             (1d + sqrt(imaginary(chi)^2+real_part(chi)^2)))
   lecc=sqrt(1d -leps^2)

   oplot,xpos[i]+r0*(1d -lecc^2)*cos(phis)/(1d +lecc*cos(phis-lpa))+$
         r0*lecc*cos(lpa),$
         ypos[i]+r0*(1d -lecc^2)*sin(phis)/(1d +lecc*cos(phis-lpa))+$
         r0*lecc*sin(lpa), thick=5d
endfor

!p.multi=0
setup_x

ps2eps,tag+'_orig.ps',/delete
ps2eps,tag+'_lens.ps',/delete

end
