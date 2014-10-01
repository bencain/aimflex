function dists, z, DA=DA, DL=DL

  omega_m=0.27d
  omega_de=0.73d
  H0=72d                        ; km/s/Mpc
  c=2.9979d5                    ; km/s

  dz=1d-3
  z=double(z)

  nz=ceil(z/dz)
  zs=z*dindgen(nz)/double(nz-1)

  D_H=c/H0

  ;Comoving transverse, luminosity and angular diameter distances
  D_com=D_H*dz*total(1d/sqrt(omega_m*(1d +zs)^3+omega_de))
  D_lum=D_com*(1d +z)
  D_ang=D_com/(1d +z)

  distances=[D_com,D_lum,D_ang]
  if keyword_set(DA) then return, D_ang
  if keyword_set(DL) then return, D_lum
  return, distances
  ; Distances are in Mpc

end

pro lensing_cosmology


  G=4.30117902d-9               ; (km/s)^2 * Mpc / M_sun
  c=2.9979d5                    ; (km/s)

  z_lens=0.9d*(dindgen(10)/9d)+0.05d
  z_src1=1.2d
  z_src2=1.0d
  z_src3=0.8d

  D_src1=dists(z_src1,/DA)
  D_src2=dists(z_src2,/DA)
  D_src3=dists(z_src3,/DA)

  D_lens=dblarr(10)
  for i=0,9 do D_lens[i]=dists(z_lens[i],/DA)

  D_ls1=D_src1 - ((1d +z_lens)/(1d +z_src1))*D_lens
  D_ls2=D_src2 - ((1d +z_lens)/(1d +z_src2))*D_lens
  D_ls3=D_src3 - ((1d +z_lens)/(1d +z_src3))*D_lens

  print,D_ls1,D_ls2,D_ls3

  ok1=where(D_ls1 gt 0d)
  ok2=where(D_ls2 gt 0d)
  ok3=where(D_ls3 gt 0d)

  sigma_crit1=(c^2/(4d*!dpi*G))*D_src1/(D_lens*D_ls1)/1d14
  sigma_crit2=(c^2/(4d*!dpi*G))*D_src2/(D_lens*D_ls2)/1d14
  sigma_crit3=(c^2/(4d*!dpi*G))*D_src3/(D_lens*D_ls3)/1d14

  file='lensing_cosmology.ps'
  t=5d

  setup_ps,file
  plot,z_lens[ok2],sigma_crit2[ok2],xtitle='Lens Redshift',$
       ytitle='Critical Surface Density (1e14 M_sun/Mpc^2)',$
       xmargin=[17,10], ymargin=[8,8],$
       charsize=2.5,charthick=t,thick=t,xthick=t,ythick=t
  oplot,z_lens[ok1],sigma_crit1[ok1],thick=t,linestyle=2
  oplot,z_lens[ok3],sigma_crit3[ok3],thick=t,linestyle=3
  xyouts,0.25,0.5,'Source redshift = 1.2, 1.0, 0.8',/norm, charsize=2.5
  setup_x
  ps2eps,file,/delete

end

