FUNCTION PENG_POWERLAW_NEW, POSITIONS, FLEXIONSCALE=FLEXIONSCALE,$
                            KAPPA=KAPPA,POLAR=POLAR

; This is a redo of the Peng lensmodel.

;Centered at
;
;Xray => (13:11:29.533,-01:20:29.67)
;WL   => (13:11:29.520,-01:20:27.59)
;diff    (00:00:00.013, 00:00:02.08)
; We take the mean for the center and convert to pixels from the HST image

pos_size=size(positions,/dimensions)
if n_elements(pos_size) gt 1 then nobj=pos_size[0] else nobj=1

if not keyword_set(flexionscale) then flexionscale=1d3

center=[1979.281d,2122.3345d] ; in pixels

if not keyword_set(polar) then begin
   x=positions[*,0]-center[0]
   y=positions[*,1]-center[1]
endif else begin
   ; Assume that the positions are in r-phi format
   if nobj gt 1 then begin
      x=positions[*,0]*cos(positions[*,1])
      y=positions[*,0]*sin(positions[*,1])
   endif else begin
      x=positions[0]*cos(positions[1])
      y=positions[0]*sin(positions[1])
   endelse
endelse

cmplx_pos=dcomplex(x,y)

lmodel_data=read_data('~/idl/peng_sigma_vs_R/kappa_combined.dat',/quiet)

; Radius is in column 0, Kappa is in column 3
rdata=lmodel_data[*,0]
kdata=lmodel_data[*,3]

; This is where the data comes from, idexed.
k_sl=lindgen(16)
k_wl=lindgen(10)+16

k_lens=lindgen(26)

k_xnpar=lindgen(9)+26
k_xpar=lindgen(36)+35

k_xray=lindgen(45)+26

; Convert radius to pixels with the scale 0.05 arcsec/pixel
rdata/=0.05d

; For now, we'll just use the X-ray data
use=k_xray

; Convert to logs
logRpos=0.5d*alog(x^2+y^2)

logK=alog(kdata[use])
logR=alog(rdata[use])

order=sort(logR)
logR=logR[order]
logK=logK[order]

hi=where(logR gt alog(100d/0.05d))
lo=where(logR lt alog(20d/0.05d))
mid=set_difference(set_difference(lindgen(n_elements(logR)),lo),hi)

; The fits give coefficients for x^0, x^1, x^2, etc in that order.
fit_lo=poly_fit(logR[lo],logK[lo],1)
fit_mid=poly_fit(logR[mid],logK[mid],1)
fit_hi=poly_fit(logR[hi],logK[hi],1)

; Note that this is the exponent for kappa as a powerlaw.  However, we
; have:
;
;  kappa = 2 * A * k^2 * (theta * conj(theta))^(k-1)
;
; so A = 2 * exp(f[0])/(f[1] + 2)^2 and k = f[1]/2 + 1

;print,transpose(fit_lo)
;print,transpose(fit_mid)
;print,transpose(fit_hi)

; Find the exact turnover points
turn1=fit_lo-fit_mid
turn1=-turn1[0]/turn1[1]

turn2=fit_hi-fit_mid
turn2=-turn2[0]/turn2[1]

;print,turn1,turn2

; Split the positions into their segments.
pos_lo=where(logRpos lt turn1,nlo)
pos_hi=where(logRpos gt turn2,nhi)
pos_mid=$
   set_difference(set_difference(lindgen(n_elements(logRpos)),pos_lo),$
                  pos_hi,count=nmid)

lmodel=dblarr(nlo+nmid+nhi,8)
lmodel[*,0:1]=positions

kappa=dblarr(nlo+nmid+nhi)

if nlo gt 0 then begin
   k=fit_lo[1]/2d + 1d
   A=exp(fit_lo[0])/(2d*k^2)
;   print,a,k
   pl_lm=powerlaw_lensmodel(cmplx_pos[pos_lo],A, k,$
                            flexionscale=flexionscale,/include_kappa)
   if nlo gt 1 then begin
      lmodel[pos_lo,2:7]=pl_lm[*,2:7] 
      kappa[pos_lo]=pl_lm[*,8]
   endif else begin
      lmodel[pos_lo,2:7]=pl_lm[2:7]
      kappa[pos_lo]=pl_lm[8]
   endelse

endif
if nmid gt 0 then begin
   k=fit_mid[1]/2d + 1d
   A=exp(fit_mid[0])/(2d*k^2)
;   print,a,k
   pl_lm=powerlaw_lensmodel(cmplx_pos[pos_mid],A,k,$
                            flexionscale=flexionscale,/include_kappa)
   if nmid gt 1 then begin
      lmodel[pos_mid,2:7]=pl_lm[*,2:7] 
      kappa[pos_mid]=pl_lm[*,8]
   endif else begin
      lmodel[pos_mid,2:7]=pl_lm[2:7]
      kappa[pos_mid]=pl_lm[8]
   endelse
endif
if nhi gt 0 then begin
   k=fit_hi[1]/2d + 1d
   A=exp(fit_hi[0])/(2d*k^2)
;   print,a,k
   pl_lm=powerlaw_lensmodel(cmplx_pos[pos_hi],A,k,$
                            flexionscale=flexionscale,/include_kappa)
   if nhi gt 1 then begin
      lmodel[pos_hi,2:7]=pl_lm[*,2:7] 
      kappa[pos_hi]=pl_lm[*,8]
   endif else begin
      lmodel[pos_hi,2:7]=pl_lm[2:7]
      kappa[pos_hi]=pl_lm[8]
   endelse

endif


if nobj eq 1 then lmodel=transpose(lmodel)

return,lmodel
end
