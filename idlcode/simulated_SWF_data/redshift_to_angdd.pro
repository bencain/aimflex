function one_over_sqrt_E,z

  ; Fixed flat concordance cosmology
  omega_matter=0.27d;0.272d
  omega_lambda=0.73d;0.728d

  return,1d/sqrt(omega_lambda + omega_matter*(1d + z)^3)
end

;;;;;;;;;;;;;;;;;;;;;;;

function redshift_to_angdd,z

; Calculate angular diameter distance between the observer and
; the object(s), in units of the Hubble distance
  DC=qromb('one_over_sqrt_E',0d*z,z,K=8,/double,JMAX=20)

  DA=DC/(1d + z)

  return,DA
end
