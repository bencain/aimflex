pro test_nie_nis


  nobj=1d6

  span=5d

  theta=span*(randomu(seed,nobj,2)-0.5d)

  ; Set up the NIE lens data
  ctr=[0d,0d]                   ; offset the NIE (in arcmin)
  sigma=1d                      ; Velocity Dispersion, in 10^3 km/s
  kappa0=0d                     ; Mass sheet offset
  theta_c=1d-4/60d              ; Core radius, in arcmin
  epsilon=1d-5                  ; Ellipticity
  position_angle=0d        ; Position angle (in radians)
  z_lens=0.2d                   ; Lens redshift

  pars=[ctr,sigma,kappa0,theta_c,epsilon,position_angle,z_lens]

  nie=nie_lens(theta, pars)

  r=sqrt(total(theta^2,2,/double))

  plot,r,sqrt(total(nie[*,8:9]^2,2)),/psym,/xlog,/ylog


;  beta=[0.0d,0.000d]
;  contour,(theta[*,0]-beta[0]-nie[*,1])^2 + (theta[*,1]-beta[1]-nie[*,1])^2, theta[*,0],theta[*,1],/irregular,nlevels=30,/downhill


end
