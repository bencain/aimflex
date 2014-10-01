function one_over_sqrt_E,z

  ; Fixed flat concordance cosmology
  omega_matter=0.27d;0.272d
  omega_lambda=0.73d;0.728d

  return,1d/sqrt(omega_lambda + omega_matter*(1d + z)^3)
end

;;;;;;;;;;;;;;;;;;;;;;;

function redshift_to_weight,z_src,z_lens
  
  z_max=1d4

  if n_elements(z_lens) gt n_elements(z_src) then $
     infty=0d*z_lens + z_max else infty=0d*z_src + z_max
  
  
; Calculate comoving line-of-sight distances between the observer and
; the source/lens and to infinity (in units of the Hubble distance).

  DL=redshift_to_angdd(z_lens)
  DS=redshift_to_angdd(z_src)
  DLS=DS - ((1d + z_lens)/(1d + z_src))*DL

  DI=redshift_to_angdd(infty)
  DLI=DI - ((1d + z_lens)/(1d + infty))*DL

  unlensed=where(z_src lt z_lens,ct)

  weight=(DLS/DS)/(DLI/DI)

  if ct gt 0 then weight[unlensed]=0d

  return,weight
end
