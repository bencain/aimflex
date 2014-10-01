function Q_FILT, X, L, R

A=(4d/sqrt(!dpi))*(gamma(3.5d + L)/gamma(3d + L))
;Q=A*((2d + L)/(2d*!dpi))*(1d - (abs(x)/abs(r))^2)^(1d + L)
Q=-A*((2d + L)/(2d*!dpi))*x*(1d - (abs(x)/abs(r))^2)^(1d + L)

return,Q
end

;;;;;;;;;;;;;;;;;;;;;;;;

function APERTURE_MASS, X, F, L_FILT, R_FILT, $
                        NX=NX, XMIN=XMIN, XMAX=XMAX, $
                        NY=NY, YMIN=YMIN, YMAX=YMAX, $
                        ABS_FE=ABS_FE, VERBOSE=VERBOSE,$
                        WEIGHTS=WEIGHTS, BMODE=BMODE

if not keyword_set(nx) then nx=1000L
nx=long(nx)
if not keyword_set(xmin) then xmin=min(real_part(x))
if not keyword_set(xmax) then xmax=max(real_part(x))
dx=(xmax-xmin)/double(nx)

if not keyword_set(ny) then ny=1000L
ny=long(ny)
if not keyword_set(ymin) then ymin=min(imaginary(x))
if not keyword_set(ymax) then ymax=max(imaginary(x))
dy=(ymax-ymin)/double(ny)

if not keyword_set(weights) then weights=dblarr(n_elements(x))+1d

mgrid=dblarr(nx,ny)
x0=dcomplex(dx*((dindgen(nx,ny) mod double(nx))+0.5d) + xmin,$
            dy*(double(lindgen(nx,ny)/long(ny))+0.5d) + ymin)

ng=nx*ny
tail=strcompress('/'+string(ng),/remove_all)
starttime=systime()

; If BMODE is set, then calculate the b-mode signal
if keyword_set(bmode) then ph=exp(dcomplex(0d,!dpi/2d)) else ph=1d

; Allow for either signed Fe or just the absolute value of the radial
; component to be used.
if keyword_set(abs_fe) then begin
   for i=0L,ng-1 do begin
      r=abs(x-x0[i])
      Fe=real_part(ph*F*conj(x-x0[i])/r) 
                                ; project all flexions to get radial
                                ; (emode) parts with respect to x0
      use=where(r lt r_filt,nuse)
      if nuse gt 0 then mgrid[i]=total(abs(Fe[use])*q_filt(r[use],l_filt,r_filt)*weights[use],/double)/$
                                 total(weights[use],/double)
   endfor
endif else begin
   for i=0L,ng-1 do begin
      r=abs(x-x0[i])
      Fe=real_part(F*conj(x-x0[i])/r) 
                                ; project all flexions to get radial
                                ; (emode) parts with respect to x0
      use=where(r lt r_filt,nuse)
      if nuse gt 0 then $
         mgrid[i]=total(Fe[use]*q_filt(r[use],l_filt,r_filt)*weights[use],/double)/$
                  total(weights[use],/double)
   endfor
endelse

if keyword_set(verbose) then begin
   print,'Start:...'+starttime
   print,'End:.....'+systime()
endif


return,mgrid

end
