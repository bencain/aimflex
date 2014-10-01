function aim_input_shear,pos,pars

; assume pars= x0,y0,theta_e,theta_c and we have a NIS lens 

shear=pos*0d

if n_elements(pars) ne 3 then return,shear

x=dcomplex(pos[*,0],pos[*,1]) - dcomplex(pars[0],pars[1])
x_e=pars[2]
x_c=pars[3]

; alpha=x_e*(sqrt(|x|^2+x_c^2)-x_c)*(x/|x|^2)

kappa=0.5d*x_e/sqrt(abs(x)^2 + x_c^2)

gamma=kappa*(2d*x_c*sqrt(abs(x)^2 + x_c^2) - abs(x)^2 -2d*x_c^2)*x^2/abs(x)^4

g=gamma/(1d - kappa)

shear[*,0]=real_part(g)
shear[*,1]=imaginary(g)

return,shear
end



  if not keyword_set(nfw_scale) then nfw_scale=1d
  if not keyword_set(nfw_ctr) then nfw_ctr=[3130d,2865d] ; Center of the lens in the image

  z_lens=0.183d
  d_ang=cosmo_dist(0.183)       ; Distance to A1689 in h_70^-1 Mpc

  M200=1.31d15*(0.7)^(-1)       ; From Peng et al. 2010
  c200=9.9d
  H0=70d                        ; km/s/Mpc
  Eofz=sqrt(0.3d*(1d + z_lens)^3 + 0.7d)
  Gconst=4.301d-9               ; (km/s)^2*Mpc/M_sun
  r200=(0.01d*Gconst*M200/(H0*Eofz)^2)^(1d/3d) ; Mpc

  r_s=r200/c200

  theta_s=(r_s/d_ang)*(180d/!dpi)*(3600d/0.05d) ; in HST pix

  a1689_lensmodel=nfw_lensmodel(cat[*,1:2],z_lens,c200,theta_s,$
                                ctr=nfw_ctr,force_scale=nfw_scale)

  if not keyword_set(detJacmin) then detJacmin=0.2d
  if not keyword_set(detJacmax) then detJacmax=10d
                                ; Limit on when to not use the shear
                                ; model and instead fix the shear to
                                ; zero.  This decision is based on the
                                ; determinant of the lensing Jacobian.


