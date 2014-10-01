pro ffttest

twodim=1

if twodim then begin
   imsize=[1001,1001]
   
   k0=2d*!dpi*dcomplex(1d/50d,0d)

   x=dcomplex((dindgen(imsize) mod imsize[0]),$
              double(lindgen(imsize)/long(imsize[0])))

   f=sin(0.5d*(k0*conj(x)+conj(k0)*x))
   fprime=k0*cos(0.5d*(k0*conj(x)+conj(k0)*x))

   kx=(2d*!dpi/imsize[0])*shift((dindgen(imsize) mod imsize[0]) - $
                                0.5d*(imsize[0]+(imsize[0] mod 2)) + 1,$
                                0.5d*(imsize[0]-(imsize[0] mod 2)) + 1,0)
   ky=(2d*!dpi/imsize[1])*shift(double(lindgen(imsize)/long(imsize[0])) - $
                                0.5d*(imsize[1]+(imsize[1] mod 2)) + 1,$
                                0,0.5d*(imsize[1]-(imsize[1] mod 2)) + 1)

   k=dcomplex(kx,ky)
endif else begin

   k0=2d*!dpi/50d
   s=101
   k=(2d*!dpi/s)*shift(dindgen(s) - $
                       0.5d*(s+(s mod 2)) + 1,$
                       0.5d*(s-(s mod 2)) + 1)

   x=dindgen(s)
   f=sin(k0*x)
   fprime=k0*cos(k0*x)
endelse

fbar=fft(f)
fprimebar=fft(fprime)
gbar=k*dcomplex(0,1)*fbar
g=fft(gbar,/inverse)

plot_hist,abs(g-fprime)


end
