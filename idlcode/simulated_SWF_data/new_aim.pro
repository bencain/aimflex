function new_aim,theta,pars,bkg=bkg

; Assume that theta is a complex position in coordinates centered on
; the image primary ray center beta(theta=0)=0


S=pars[0]
a=pars[1]
n=pars[2]

shear=dcomplex(pars[3],pars[4])
fflex=dcomplex(pars[5],pars[6])
gflex=dcomplex(pars[7],pars[8])


b_wl=theta-shear*conj(theta)
b_fl=-0.25d*conj(fflex)*theta^2 - 0.5d*fflex*theta*conj(theta) - 0.25d*gflex*conj(theta)^2

I=S*exp(-(abs(b_wl)/a)^(1d/n))*(1d - (abs(b_wl)/a)^(1d/n)*(b_wl*conj(b_fl)+conj(b_wl)*b_fl)/(2d*abs(b_wl)^2))


return,real_part(I)

end
