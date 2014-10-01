pro final_out

  ra=104.5665559d               ; DD
  dec=-55.9434082d              ; DD
  hdrfile='~/idl/pro/simulated_SWF_data/1e0657R_avg_wcs_hdr.fits'

  density=150d                  ; arcmin^-2
  field=3.5d                    ; arcmin

; Set up the NIE lens data
  ctr=[0d,0d]                   ; offset the NIE (in arcmin)
  sigma=1d                      ; Velocity Dispersion, in 10^3 km/s
  kappa0=0d                     ; Mass sheet offset
  theta_c=4d-2                  ; Core radius, in arcmin
  q=0.75d
  epsilon=(1d - q)/(1d + q)     ; Ellipticity
  position_angle=0d             ; Position angle (in radians)
  z_lens=0.25d                  ; Lens redshift
  DLI_DI = ( redshift_to_angdd(1d4) - (1d + z_lens)*redshift_to_angdd(z_lens)/(1d + 1d4) )/redshift_to_angdd(1d4)
                                ; Lensing weight

;  nie1=[ctr[0],ctr[1], $
;        sigma,kappa0,theta_c, $
;        epsilon,position_angle,DLI_DI] ;horizontal ellipse
  
;  nie2=nie1                     ;angled ellipse
;  nie2[6]=0.2
  
;  nis=nie1                      ;Sphere
;  nis[5]=0d
;  nis[6]=0d

;Substructures
;  sub1=[0.15,0.0, 0.75d, 0d, 5d-3, 0d,0d, DLI_DI]
;  sub1=[0.2,0.0, 0.5d, 0d, 5d-3, 0d,0d, DLI_DI]
;  sub3=[0.3,0.0, 0.5d, 0d, 5d-3, 0d,0d, DLI_DI]
;  sub4=[0.4,0.0, 0.5d, 0d, 5d-3, 0d,0d, DLI_DI]
  

;;;;;;;;;;;;;;;;;;;;;

;  nie=[ctr,$
;       1.00000000E+00,  0.00000000E+00,  4.00000000E-02,  0.00000000E+00,  0.00000000E+00, $ ; from the lensmodel file...
;       DLI_DI]


  minx=-(field/2d)
  miny=-(field/2d)


  dir='~/analysis_code/swunited_rewrite/simdata/perf2/'
  
  true=read_data(dir+'perf2.dat')

;;;;

  readcol,dir+'GR_interp_final.dat',x,y,k,ga1,ga2,F1,F2,G1,G2
  ngr=n_elements(x)
  int_gr=dblarr(ngr,9)
  int_gr[*,0]=x + minx
  int_gr[*,1]=y + miny

  int_gr[*,2]=k
  int_gr[*,3]=ga1
  int_gr[*,4]=ga2

  int_gr[*,5]=F1
  int_gr[*,6]=F2
  int_gr[*,7]=G1
  int_gr[*,8]=G2

  theta_gr=dblarr(n_elements(x),2)
  theta_gr[*,0]=x + minx
  theta_gr[*,1]=y + miny
  delvarx,x,y,k,ga1,ga2,F1,F2,G1,G2

;;;;

  readcol,dir+'FL_interp_final.dat',x,y,k,ga1,ga2,F1,F2,G1,G2
  nfl=n_elements(x)
  int_fl=dblarr(nfl,9)
  int_fl[*,0]=x + minx
  int_fl[*,1]=y + miny

  int_fl[*,2]=k
  int_fl[*,3]=ga1
  int_fl[*,4]=ga2

  int_fl[*,5]=F1
  int_fl[*,6]=F2
  int_fl[*,7]=G1
  int_fl[*,8]=G2

  theta_fl=dblarr(n_elements(x),2)
  theta_fl[*,0]=x + minx
  theta_fl[*,1]=y + miny
  delvarx,x,y,k,ga1,ga2,F1,F2,G1,G2

;;;;

  readcol,dir+'SL_interp_final.dat',x,y,a1,a2,k,ga1,ga2,det,zed,zwt
  nsl=n_elements(x)
  int_sl=dblarr(nsl,10)
  int_sl[*,0]=x; + minx
  int_sl[*,1]=y; + miny

  int_sl[*,2]=a1
  int_sl[*,3]=a2

  int_sl[*,4]=k
  int_sl[*,5]=ga1
  int_sl[*,6]=ga2

  int_sl[*,7]=det
  int_sl[*,8]=zed
  int_sl[*,9]=zwt

  theta_sl=dblarr(n_elements(x),2)
  theta_sl[*,0]=x; + minx
  theta_sl[*,1]=y; + miny
  delvarx,x,y,a1,a2,k,ga1,ga2,det,zwt,zed


  readcol,dir+'WL_interp_final.dat',x,y,k,ga1,ga2
  nwl=n_elements(x)
  int_wl=dblarr(nwl,5)
  int_wl[*,0]=x + minx
  int_wl[*,1]=y + miny

  int_wl[*,2]=k
  int_wl[*,3]=ga1
  int_wl[*,4]=ga2

  theta_wl=dblarr(n_elements(x),2)
  theta_wl[*,0]=x + minx
  theta_wl[*,1]=y + miny
  delvarx,x,y,k,ga1,ga2


  kout=dblarr(ngr+nfl+nsl+nwl,3)
  kout[*,0]=[int_gr[*,0],int_fl[*,0],int_sl[*,0],int_wl[*,0]]
  kout[*,1]=[int_gr[*,1],int_fl[*,1],int_sl[*,1],int_wl[*,1]]
  kout[*,2]=[int_gr[*,2],int_fl[*,2],int_sl[*,4],int_wl[*,2]]


  lvls=dindgen(31)/3d
  contour,kout[*,2],kout[*,0],kout[*,1],/irregular,/isotropic,levels=lvls

  slobj=where(true[*,21] gt 1)
  plot,true[slobj,2],true[slobj,3],/psym,/isotropic
;  plot,true[slobj,8],true[slobj,9],psym=4,/isotropic

end
