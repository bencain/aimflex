pro aim_par_corr

; Test the correlations between various aim parameter
; formulations.  We'll call them the physical and abstract
; parameterizations. 

  npar=1000L

; First, parameters in the physical parameterization

  alpha=randomu(seed,npar,/double)
  eps=dcomplex(randomu(seed,npar,/double)-0.5d,randomu(seed,npar,/double)-0.5d)
  
  shear=dcomplex(randomu(seed,npar,/double)-0.5d,randomu(seed,npar,/double)-0.5d)*3d
  fflex=dcomplex(randomu(seed,npar,/double)-0.5d,randomu(seed,npar,/double)-0.5d)*5d
  gflex=dcomplex(randomu(seed,npar,/double)-0.5d,randomu(seed,npar,/double)-0.5d)*5d

; Now create parameters in the abstract parameterization

  m=1+eps*conj(shear)

  a=alpha/abs(m)

  s=(eps+shear)/m
  t=(fflex-eps*gflex)/conj(m)
  u=(fflex-eps*conj(fflex))/m
  v=(gflex-eps*fflex)/m

  data=[transpose(a),transpose(alpha),$
        transpose(real_part(s)),transpose(imaginary(s)),transpose(real_part(eps)),transpose(imaginary(eps)),$
        transpose(real_part(t)),transpose(imaginary(t)),transpose(real_part(shear)),transpose(imaginary(shear)),$
        transpose(real_part(u)),transpose(imaginary(u)),transpose(real_part(fflex)),transpose(imaginary(fflex)),$
        transpose(real_part(v)),transpose(imaginary(v)),transpose(real_part(gflex)),transpose(imaginary(gflex))]
  

  cov=correlate(data)

  print, floor(2*cov)

  !p.multi=[0,3,3]
  plot,real_part(eps),imaginary(eps),/psym,xrange=[-1,1],yrange=[-1,1],charsize=2
  plot,real_part(shear),imaginary(shear),/psym,xrange=[-3,3],yrange=[-3,3],charsize=2
  plot,real_part(s),imaginary(s),/psym,xrange=[-3,3],yrange=[-3,3],charsize=2

  plot,real_part(fflex),imaginary(fflex),/psym,xrange=[-5,5],yrange=[-5,5],charsize=2
  plot,real_part(t),imaginary(t),/psym,xrange=[-5,5],yrange=[-5,5],charsize=2
  plot,real_part(u),imaginary(u),/psym,xrange=[-5,5],yrange=[-5,5],charsize=2


  plot,real_part(gflex),imaginary(gflex),/psym,xrange=[-5,5],yrange=[-5,5],charsize=2
  plot,real_part(v),imaginary(v),/psym,xrange=[-5,5],yrange=[-5,5],charsize=2
  
  plot,alpha,a,/psym,charsize=2
  !p.multi=0

end
