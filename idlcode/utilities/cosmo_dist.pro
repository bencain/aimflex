function cosmo_dist,redshift,omega_m=omega_m, omega_de=omega_de,$
                        omega_r=omega_r, dz=dz, h70=h70, eps=eps, verb=verb


; This function estimates the angular diameter distance between two
; redshifts. REDSHIFT is assumed to be a single redshift or an array
; of redshifts.  In the former case, the integration runs from zero to
; redshift.  In the latter case, the integration is from the minimum
; redshift to the maximum.  The result is given in Mpc.

; The standard OMEGA_M=0.3, OMEGA_DE=0.7 cosmology is used unless
; otherwise specified.  DZ is the redshift resolution for the
; integration.  This largely follows from Hogg (2000). 

if not keyword_set(dz) then dz=0.0001d

z2=double(max(redshift))

if not keyword_set(omega_m) then omega_m=0.2999d
if not keyword_set(omega_de) then omega_de=0.7d
if not keyword_set(omega_r) then omega_r=(1d)/10000d
if not keyword_set(h70) then h70=1d

if not keyword_set(eps) then eps=double(1e-10)

d_hubble=((2.99792458d)/(7d))*(10d)^4/h70 ; In Mpc

omega_k=(1d)-(omega_m+omega_de+omega_r)

nz2=long(z2/dz)

zs2=(z2/double(nz2-1))*dindgen(nz2)

ez_inv2=(1d)/sqrt(omega_r*(1+zs2)^4+omega_m*(1+zs2)^3+$
                  omega_k*(1+zs2)^2+omega_de)

dc2=tsum(zs2,ez_inv2)

if abs(omega_k) lt eps then begin
   D2=dc2
endif else if omega_k gt 0 then begin
   D2=sinh(dc2*sqrt(omega_k))/sqrt(omega_k)
endif else begin
   D2=sin(dc2*sqrt(abs(omega_k)))/sqrt(abs(omega_k))
endelse
; Convert to angular diameter distance and Mpc
D2/=z2+1d
D2*=d_hubble

if n_elements(redshift) gt 1 then begin
   z1=double(min(redshift))
   D1=cosmo_dist(z1, omega_m=omega_m, omega_de=omega_de,$
                 omega_r=omega_r, dz=dz, h70=h70, eps=eps) 
endif else begin
   z1=0d
   D1=0d
endelse

D12=D2*sqrt((1d)+omega_k*(((1d)+z1)*D1/d_hubble)^2)-$
    (((1d)+z1)/((1d)+z2))*D1*sqrt((1d)+omega_k*(((1d)+z2)*D2/d_hubble)^2)

c2_4piG=(1.66281534d)*(10d)^18d ; M_solar/Mpc

sigcrit=c2_4piG*D2/(max([D1,1d-10],ind)*D12)

scale1=D1*!dpi/648d ; 1 arcsec at z_1
scale2=D2*!dpi/648d ; 1 arcsec at z_2

if keyword_set(verb) then begin
   print,'*** ***'
   print,'z_1 = '+strcompress(string(z1),/remove_all)
   print,'z_2 = '+strcompress(string(z2),/remove_all)
   print,'***Cosmology***'
   print,'O_m = '+strcompress(string(omega_m),/remove_all)+', '+$
         'O_de = '+strcompress(string(omega_de),/remove_all)
   print,'O_r = '+strcompress(string(omega_r),/remove_all)+', '+$
         'O_k = '+strcompress(string(omega_k),/remove_all)
   print,'h_70 = '+strcompress(string(h70),/remove_all)
   print,'*** ***'
   print,'D_1 = '+strcompress(string(d1),/remove_all)+' h_70^-1 Mpc'
   print,'D_2 = '+strcompress(string(d2),/remove_all)+' h_70^-1 Mpc'
   print,'D_12 = '+strcompress(string(d12),/remove_all)+' h_70^-1 Mpc'
   if ind eq 1 then  print,'Sigma_crit = Infinity' else $
      print,'Sigma_crit = '+$
            strcompress(string(sigcrit),/remove_all)+' h_70 M_solar/Mpc^2'
   print,'1 arcsec at z_1 = '+strcompress(scale1,/remove_all)+' h_70^-1 kpc'
   print,'1 arcsec at z_2 = '+strcompress(scale2,/remove_all)+' h_70^-1 kpc'
   print,'*** ***'
endif

return,d12

end
