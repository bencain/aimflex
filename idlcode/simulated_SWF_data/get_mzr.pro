function get_mzr,n,seed=seed,fixedm=fixedm,fixedz=fixedz,fixedr=fixedr

; n is the number of galaxies to draw

  if not keyword_set(n) then n=1L
  if n lt 1L then n=1L

  out=dblarr(n,3)


  if not keyword_set(fixedz) then begin
; Draw the redshifts from a gamma distribution
     z0=1d/3d                   ;From Schrabback+2010
     zmin=0.3d
     zmax=5d

     nfilled=0
     while nfilled lt n do begin
        z=z0*randomn(seed,/double,gamma=3)
        if (z lt zmax) and (z gt zmin) then begin
           out[nfilled,1]=z
           nfilled++
        endif
     endwhile

;;;;;;;;;;;;;;;;
;  OVERRIDING...uniform dist from z=1 to z=3
     out[*,1]=2d*randomu(seed,n) + 1d

;;;;;;;;;;;;;;;;

  endif else out[*,1]=fixedz[0]


; For now we'll do something rough (and probably wrong) just to
; move along.  We'll assume a roughly Schechter-ish function
; for the apparent magnitude at z=1, draw from that and then adjust
; the luminosity based on the actual redshift.
  if not keyword_set(fixedm) then begin
     mmax=29
     mmin=19
     dm=(mmax-mmin)/1d4
   
;  dh=(2.9979d5/70d)             ;Hubble dist in Mpc
     dldh_z_1=redshift_to_angdd(1d)*4d ; Luminosity distant to z=1 in Hubble distances

     alpha=-1.05
     mstar=23d
     norm=0.16039766d
  
     nfilled=0
     while nfilled lt n do begin
        m=(mmax-mmin)*randomu(seed,/double) + mmin
        m_rnd=ceil(m/dm,/l64)*dm
        p=randomu(seed)
      
        p_of_m=norm*(10d)^(-0.4d*(m_rnd-mstar)*(alpha+1d))*exp(-(10d)^(-0.4d*(m_rnd-mstar)))*dm
        if p lt p_of_m then begin
           z=out[nfilled,1]
           out[nfilled,0]=m-5d*alog10(dldh_z_1/(redshift_to_angdd(z)*(1d + z)^2))
           nfilled++
        endif
     endwhile
  endif else out[*,0]=fixedm[0]

; For calculating norm
;  ms=mmax-dm*dindgen(1d4)
;  tot=total(norm*(10d)^(-0.4d*(ms-mstar)*(alpha+1d))*exp(-(10d)^(-0.4d*(ms-mstar)))*dm)
;  print,1d/tot
;  plot,ms,total(norm*(10d)^(-0.4d*(ms-mstar)*(alpha+1d))*exp(-(10d)^(-0.4d*(ms-mstar)))*dm,/cum)
; 1/tot becomes norm

  if not keyword_set(fixedr) then begin
     r0=0.5d                    ; size at z=1
     out[*,2]=r0*redshift_to_angdd(1d)/redshift_to_angdd(out[*,1])
  endif else out[*,2]=fixedr[0]

  return,out
end
