pro dawson,ntry
; Simulation for the Dawson et al. 2012 Chandra+HST proposal
  if ntry lt 10L then ntry=10L else ntry=long(ntry)

; All masses will be in 10^14 solar masses, NFW sphere profiles.

  z=0.53d
  DL=(299792.458d/70d)*redshift_to_angdd(z) ; In Mpc
  
  mass_south_dm=3.09d
  mass_south_gas=0.1d
  mass_cent_gas=0.1d
  mass_north=1.69d              ; all in 1e14 Msun

  rs_south_dm=(0.39d/DL)*(180d*3600d/!dpi) 
  rs_south_gas=(0.09d/DL)*(180d*3600d/!dpi)
  rs_cent_gas=(0.09d/DL)*(180d*3600d/!dpi)
  rs_north=(0.30d/DL)*(180d*3600d/!dpi) ; all in arcsec

  dc_south_dm=2776d
  dc_south_gas=5690d
  dc_cent_gas=5690d
  dc_north=3132d                ; all unitless


; We'll assume everything is aligned along the x-axis.  These
; are the "null" positions.  The south and north DM centers are
; assumed to be at the galaxy positions.  I'll be varying the
; south DM position about this center.

  ctr_south_dm=[-(191d/2d),0d]
  ctr_south_gas=[40d - (191d/2d),0d]
  ctr_cent_gas=[106d - (191d/2d),0d]
  ctr_north=[(191d/2d),0d]

; Make the base lens arrays
  lens_south_dm=[ctr_south_dm,$
                 rs_south_dm,$
                 dc_south_dm,z]
  lens_south_gas=[ctr_south_gas,$
                  rs_south_gas,$
                  dc_south_gas,z]
  lens_cent_gas=[ctr_cent_gas,$
                 rs_cent_gas,$
                 dc_cent_gas,z]
  lens_north=[ctr_north,$
              rs_north,$
              dc_north,z]


; Set up the flexion measurement positions in an n by n grid

  L=600d                        ; Side length, in arcsec
  n=600L                        ; pixels per side (square grid)
  
  theta=dblarr(n^2,2)
  theta[*,0]=(dindgen(n^2) mod n)/double(n) - 0.5d
  theta[*,1]=double(lindgen(n^2)/n)/double(n) - 0.5d  
  theta*=L

  theta+=!dpi/100d

; Calculate the fields from the fixed components

  fields_fixed= nfw_lens(theta,lens_south_gas,/pars_delta_c)
  fields_fixed+=nfw_lens(theta,lens_cent_gas,/pars_delta_c)
  fields_fixed+=nfw_lens(theta,lens_north,/pars_delta_c)

; The null hypothesis fields
  fields_null=fields_fixed+nfw_lens(theta,lens_south_dm,/pars_delta_c)

  aF1=dblarr(n,n)
  aF1[*]=0.5d*fields_null[*,6]/(1d - fields_null[*,3])
  aF2=dblarr(n,n)
  aF2[*]=0.5d*fields_null[*,7]/(1d - fields_null[*,3])

  mwrfits,sqrt(aF1^2 + aF2^2),'nullaF.fits'
  
; Set up output space
  aF1_out=dblarr(ntry,n^2)
  aF2_out=dblarr(ntry,n^2)
  aF_out=dblarr(ntry,n^2)
  r=dblarr(ntry)

  for i=0,ntry-1 do begin

     lens=lens_south_dm
     dx=randomn(seed,2,/double)*20d ; Vary position in a 20 arcsec gaussian

     lens[0:1]+=dx
     r[i]=sqrt(total(dx^2))

     fields=fields_fixed+nfw_lens(theta,lens,/pars_delta_c)
  
 ; Assume a 0.5 arcsec image to make dimensionless
     aF1_out[i,*]=0.5d*(fields[*,6]/(1d - fields[*,3]))
     aF2_out[i,*]=0.5d*(fields[*,7]/(1d - fields[*,3]))
  endfor

  aF_out=sqrt(aF1_out^2+aF2_out^2)

  sig_aF=dblarr(n,n)
  for i=0,n^2-1 do sig_aF[i]=stddev(aF_out[*,i])
  mwrfits,sig_aF,'sig_aF.fits'

end
