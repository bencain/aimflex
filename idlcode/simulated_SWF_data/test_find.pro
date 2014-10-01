pro test_find

  ra=104.5665559d               ; DD
  dec=-55.9434082d              ; DD
  hdrfile='pro/simulated_SWF_data/1e0657R_avg_wcs_hdr.fits'

  density=200d                  ; arcmin^-2
  field=3.5d                    ; arcmin
  wn=0                          ; noise? NOTE - I've hand set the SL and FL noise to zero, no matter what.

; Set up the NIE lens data
  ctr=[0d,0d]                   ; offset the NIE (in arcmin)
  sigma=1d                      ; Velocity Dispersion, in 10^3 km/s
  kappa0=0d                     ; Mass sheet offset
  theta_c=5d-3                  ; Core radius, in arcmin
  q=0.5d
  epsilon=(1d - q)/(1d + q)     ; Ellipticity
  position_angle=0d        ; Position angle (in radians)
  z_lens=0.45d                  ; Lens redshift
  DLI_DI = ( redshift_to_angdd(1d4) - (1d + z_lens)*redshift_to_angdd(z_lens)/(1d + 1d4) )/redshift_to_angdd(1d4)

  nie0=[ctr,sigma,kappa0,theta_c,epsilon,position_angle,DLI_DI]
  nis0=[ctr,sigma,kappa0,theta_c,DLI_DI]

  Z=0.7d
  err=1d-10
  span=5d
  
  theta=1d*(randomu(seed,5000,2,/double)-0.5d)  
  fields=my_lens(theta,nie0,z_weight=Z)

  for i=0,99 do begin

     beta=0.05d*(randomu(seed,2,/double)-0.5d)

     imgs=find_mult_images_triangles(beta,err,Z,nie0,span=5d,ngrid=20,nimages=nimgs)

     if nimgs gt 1 then begin
        contour,fields[*,3],theta[*,0],theta[*,1],/isotropic,/irregular
        oplot,imgs[*,0],imgs[*,1],/psym,color=fsc_color('red')
        oplot,[beta[0]],[beta[1]],/psym,color=fsc_color('blue')
        read,blar
     endif
  endfor
  

end
