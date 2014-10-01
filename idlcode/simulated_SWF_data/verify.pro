pro verify
;;;;;;;;;;;;;;;;;;;;;;


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

  nie1=[ctr[0],ctr[1], $
        sigma,kappa0,theta_c, $
        epsilon,position_angle,DLI_DI] ;horizontal ellipse
  
  nie2=nie1                     ;angled ellipse
  nie2[6]=0.2
  
  nis=nie1                      ;Sphere
  nis[5]=0d
  nis[6]=0d

;Substructures
  sub1=[0.15,0.0, 0.75d, 0d, 5d-3, 0d,0d, DLI_DI]
  sub1=[0.2,0.0, 0.5d, 0d, 5d-3, 0d,0d, DLI_DI]
  sub3=[0.3,0.0, 0.5d, 0d, 5d-3, 0d,0d, DLI_DI]
  sub4=[0.4,0.0, 0.5d, 0d, 5d-3, 0d,0d, DLI_DI]
  

;;;;;;;;;;;;;;;;;;;;;

  nie=[ctr,$
       1.00000000E+00,  0.00000000E+00,  4.00000000E-02,  0.00000000E+00,  0.00000000E+00, $ ; from the lensmodel file...
       DLI_DI]


  minx=-(field/2d)
  miny=-(field/2d)

  dir='~/analysis_code/swunited_rewrite/simdata/perf2/'
  data=read_data(dir+'perf.dat')

  readcol,dir+'GR_interp.dat',x,y,k,ga1,ga2,F1,F2,G1,G2
 
  int_gr=dblarr(n_elements(x),9)
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

;;;;;;;;;;;;;;; flexion excluded for now...
if 0 then begin
  readcol,dir+'FL_interp.dat',x,y,k,ga1,ga2,F1,F2,G1,G2
  int_fl=dblarr(n_elements(x),9)
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
endif

  readcol,dir+'SL_interp.dat',x,y,a1,a2,k,ga1,ga2,det,zed,zwt
  int_sl=dblarr(n_elements(x),10)
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


  readcol,dir+'WL_interp.dat',x,y,k,ga1,ga2
  int_wl=dblarr(n_elements(x),5)
  int_wl[*,0]=x + minx
  int_wl[*,1]=y + miny

  int_wl[*,2]=k
  int_wl[*,3]=ga1
  int_wl[*,4]=ga2

  theta_wl=dblarr(n_elements(x),2)
  theta_wl[*,0]=x + minx
  theta_wl[*,1]=y + miny
  delvarx,x,y,k,ga1,ga2

;;;;;;;;;;;;;;;;;;;;;;

  model_gr=nie_lens(theta_gr,nie)
  model_wl=nie_lens(theta_wl,nie)
;  model_fl=nie_lens(theta_fl,nie)
  model_sl=nie_lens(theta_sl,nie)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Grid plots

  set_plot,'PS'
  !p.font=0
  device,filename='verify_GRID_kappa.eps',$
         xsize=6,ysize=6,/inches,/decomposed,/color,bits_per_pixel=8,/encapsulated
  plot, model_gr[*,3],int_gr[*,2],/psym,$
        title='kappa grid',xtitle='model value',ytitle='interpreted value'
  oplot,[-500,500],[-500,500]
  device,/close
  set_plot,'X'

  set_plot,'PS'
  !p.font=0
  device,filename='verify_GRID_gamma1.eps',$
         xsize=6,ysize=6,/inches,/decomposed,/color,bits_per_pixel=8,/encapsulated
  plot, model_gr[*,4],int_gr[*,3],/psym,$
        title='gamma1 grid',xtitle='model value',ytitle='interpreted value'
  oplot,[-500,500],[-500,500]
  device,/close
  set_plot,'X'

  set_plot,'PS'
  !p.font=0
  device,filename='verify_GRID_gamma2.eps',$
         xsize=6,ysize=6,/inches,/decomposed,/color,bits_per_pixel=8,/encapsulated
  plot, model_gr[*,5],int_gr[*,4],/psym,$
        title='gamma2 grid',xtitle='model value',ytitle='interpreted value'
  oplot,[-500,500],[-500,500]
  device,/close
  set_plot,'X'

  set_plot,'PS'
  !p.font=0
  device,filename='verify_GRID_mu.eps',$
         xsize=6,ysize=6,/inches,/decomposed,/color,bits_per_pixel=8,/encapsulated

  mu_gr_m = 1d/((1d - model_gr[*,3])^2 - (model_gr[*,4]^2 + model_gr[*,5]^2))
  mu_gr_i = 1d/((1d - int_gr[*,2])^2 - (int_gr[*,3]^2 + int_gr[*,4]^2))

  plot, mu_gr_m,mu_gr_i,/psym,$
        title='mu grid',xtitle='model value',ytitle='interpreted value'
  oplot,[-500,500],[-500,500]
  device,/close
  set_plot,'X'

  set_plot,'PS'
  !p.font=0
  device,filename='verify_GRID_fflex1.eps',$
         xsize=6,ysize=6,/inches,/decomposed,/color,bits_per_pixel=8,/encapsulated
  plot, model_gr[*,6],int_gr[*,5],/psym,$
        title='fflex1 grid',xtitle='model value',ytitle='interpreted value'
  oplot,[-500,500],[-500,500]
  device,/close
  set_plot,'X'

  set_plot,'PS'
  !p.font=0
  device,filename='verify_GRID_fflex1.eps',$
         xsize=6,ysize=6,/inches,/decomposed,/color,bits_per_pixel=8,/encapsulated
  plot, model_gr[*,6],int_gr[*,5],/psym,$
        title='fflex1 grid',xtitle='model value',ytitle='interpreted value'
  oplot,[-500,500],[-500,500]
  device,/close
  set_plot,'X'

  set_plot,'PS'
  !p.font=0
  device,filename='verify_GRID_fflex2.eps',$
         xsize=6,ysize=6,/inches,/decomposed,/color,bits_per_pixel=8,/encapsulated
  plot, model_gr[*,7],int_gr[*,6],/psym,$
        title='fflex2 grid',xtitle='model value',ytitle='interpreted value'
  oplot,[-500,500],[-500,500]
  device,/close
  set_plot,'X'

  set_plot,'PS'
  !p.font=0
  device,filename='verify_GRID_gflex1.eps',$
         xsize=6,ysize=6,/inches,/decomposed,/color,bits_per_pixel=8,/encapsulated
  plot, model_gr[*,8],int_gr[*,7],/psym,$
        title='gflex1 grid',xtitle='model value',ytitle='interpreted value'
  oplot,[-500,500],[-500,500]
  device,/close
  set_plot,'X'

  set_plot,'PS'
  !p.font=0
  device,filename='verify_GRID_gflex2.eps',$
         xsize=6,ysize=6,/inches,/decomposed,/color,bits_per_pixel=8,/encapsulated
  plot, model_gr[*,9],int_gr[*,8],/psym,$
        title='gflex2 grid',xtitle='model value',ytitle='interpreted value'
  oplot,[-500,500],[-500,500]
  device,/close
  set_plot,'X'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; SL point plots

  set_plot,'PS'
  !p.font=0
  device,filename='verify_SL_alpha1.eps',$
         xsize=6,ysize=6,/inches,/decomposed,/color,bits_per_pixel=8,/encapsulated
  plot, model_sl[*,1],int_sl[*,2],/psym,$
        title='alpha1 SL',xtitle='model value',ytitle='interpreted value'
  oplot,[-500,500],[-500,500]
  device,/close
  set_plot,'X'

  set_plot,'PS'
  !p.font=0
  device,filename='verify_SL_alpha2.eps',$
         xsize=6,ysize=6,/inches,/decomposed,/color,bits_per_pixel=8,/encapsulated
  plot, model_sl[*,2],int_sl[*,3],/psym,$
        title='alpha2 SL',xtitle='model value',ytitle='interpreted value'
  oplot,[-500,500],[-500,500]
  device,/close
  set_plot,'X'

  set_plot,'PS'
  !p.font=0
  device,filename='verify_SL_kappa.eps',$
         xsize=6,ysize=6,/inches,/decomposed,/color,bits_per_pixel=8,/encapsulated
  plot, model_sl[*,3],int_sl[*,4],/psym,$
        title='kappa SL',xtitle='model value',ytitle='interpreted value'
  oplot,[-500,500],[-500,500]
  device,/close
  set_plot,'X'

  set_plot,'PS'
  !p.font=0
  device,filename='verify_SL_gamma1.eps',$
         xsize=6,ysize=6,/inches,/decomposed,/color,bits_per_pixel=8,/encapsulated
  plot, model_sl[*,4],int_sl[*,5],/psym,$
        title='gamma1 SL',xtitle='model value',ytitle='interpreted value'
  oplot,[-500,500],[-500,500]
  device,/close
  set_plot,'X'

  set_plot,'PS'
  !p.font=0
  device,filename='verify_SL_gamma2.eps',$
         xsize=6,ysize=6,/inches,/decomposed,/color,bits_per_pixel=8,/encapsulated
  plot, model_sl[*,5],int_sl[*,6],/psym,$
        title='gamma2 SL',xtitle='model value',ytitle='interpreted value'
  oplot,[-500,500],[-500,500]
  device,/close
  set_plot,'X'

  Zwt_m=redshift_to_weight(int_sl[*,8],z_lens)

  mu_sl_m_field = 1d/((1d - model_sl[*,3])^2 - (model_sl[*,4]^2 + model_sl[*,5]^2))
  mu_sl_i_field = 1d/((1d - int_sl[*,4])^2 - (int_sl[*,5]^2 + int_sl[*,6]^2))

  mu_sl_m_det = 1d/((1d - Zwt_m*model_sl[*,3])^2 - Zwt_m^2*(model_sl[*,4]^2 + model_sl[*,5]^2))
  mu_sl_i_det=1d/int_sl[*,7]

  set_plot,'PS'
  !p.font=0
  device,filename='verify_SL_mu_kg.eps',$
         xsize=6,ysize=6,/inches,/decomposed,/color,bits_per_pixel=8,/encapsulated
  plot, mu_sl_m_field,mu_sl_i_field,/psym,$
        title='mu SL inf src',xtitle='model value',ytitle='interpreted value'
  oplot,[-500,500],[-500,500]
  device,/close
  set_plot,'X'

  set_plot,'PS'
  !p.font=0
  device,filename='verify_SL_mu_det.eps',$
         xsize=6,ysize=6,/inches,/decomposed,/color,bits_per_pixel=8,/encapsulated
  plot, mu_sl_m_det,mu_sl_i_det,/psym,$
        title='mu SL finite src z',xtitle='model value',ytitle='interpreted value'
  oplot,[-500,500],[-500,500]
  device,/close
  set_plot,'X'


  set_plot,'PS'
  !p.font=0
  device,filename='verify_SL_zwt.eps',$
         xsize=6,ysize=6,/inches,/decomposed,/color,bits_per_pixel=8,/encapsulated
  plot, Zwt_m,int_sl[*,9],/psym,$
        title='Z(z)',xtitle='model value',ytitle='interpreted value'
  oplot,[-500,500],[-500,500]
  device,/close
  set_plot,'X'


  nsys=n_elements(int_sl[*,0])/4L
  chisq=0d
  sig=0.9d/60d
     


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; WL plots

  set_plot,'PS'
  !p.font=0
  device,filename='verify_WL_kappa.eps',$
         xsize=6,ysize=6,/inches,/decomposed,/color,bits_per_pixel=8,/encapsulated
  plot, model_wl[*,3],int_wl[*,2],/psym,$
        title='kappa WL',xtitle='model value',ytitle='interpreted value'
  oplot,[-500,500],[-500,500]
  device,/close
  set_plot,'X'

  set_plot,'PS'
  !p.font=0
  device,filename='verify_WL_gamma1.eps',$
         xsize=6,ysize=6,/inches,/decomposed,/color,bits_per_pixel=8,/encapsulated
  plot, model_wl[*,4],int_wl[*,3],/psym,$
        title='gamma1 WL',xtitle='model value',ytitle='interpreted value'
  oplot,[-500,500],[-500,500]
  device,/close
  set_plot,'X'

  set_plot,'PS'
  !p.font=0
  device,filename='verify_WL_gamma2.eps',$
         xsize=6,ysize=6,/inches,/decomposed,/color,bits_per_pixel=8,/encapsulated
  plot, model_wl[*,5],int_wl[*,4],/psym,$
        title='gamma2 WL',xtitle='model value',ytitle='interpreted value'
  oplot,[-500,500],[-500,500]
  device,/close
  set_plot,'X'

  set_plot,'PS'
  !p.font=0
  device,filename='verify_WL_mu.eps',$
         xsize=6,ysize=6,/inches,/decomposed,/color,bits_per_pixel=8,/encapsulated

  mu_wl_m = 1d/((1d - model_wl[*,3])^2 - (model_wl[*,4]^2 + model_wl[*,5]^2))
  mu_wl_i = 1d/((1d - int_wl[*,2])^2 - (int_wl[*,3]^2 + int_wl[*,4]^2))

  plot, mu_wl_m,mu_wl_i,/psym,$
        title='mu WL',xtitle='model value',ytitle='interpreted value'
  oplot,[-500,500],[-500,500]
  device,/close
  set_plot,'X'


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; FL plots
;;;;;;;;;;;;;;;
if 0 then begin

  set_plot,'PS'
  !p.font=0
  device,filename='verify_FL_kappa.eps',$
         xsize=6,ysize=6,/inches,/decomposed,/color,bits_per_pixel=8,/encapsulated
  plot, model_fl[*,3],int_fl[*,2],/psym,$
        title='kappa FL',xtitle='model value',ytitle='interpreted value'
  oplot,[-500,500],[-500,500]
  device,/close
  set_plot,'X'

  set_plot,'PS'
  !p.font=0
  device,filename='verify_FL_gamma1.eps',$
         xsize=6,ysize=6,/inches,/decomposed,/color,bits_per_pixel=8,/encapsulated
  plot, model_fl[*,4],int_fl[*,3],/psym,$
        title='gamma1 FL',xtitle='model value',ytitle='interpreted value'
  oplot,[-500,500],[-500,500]
  device,/close
  set_plot,'X'

  set_plot,'PS'
  !p.font=0
  device,filename='verify_FL_gamma2.eps',$
         xsize=6,ysize=6,/inches,/decomposed,/color,bits_per_pixel=8,/encapsulated
  plot, model_fl[*,5],int_fl[*,4],/psym,$
        title='gamma2 FL',xtitle='model value',ytitle='interpreted value'
  oplot,[-500,500],[-500,500]
  device,/close
  set_plot,'X'

  set_plot,'PS'
  !p.font=0
  device,filename='verify_FL_mu.eps',$
         xsize=6,ysize=6,/inches,/decomposed,/color,bits_per_pixel=8,/encapsulated

  mu_fl_m = 1d/((1d - model_fl[*,3])^2 - (model_fl[*,4]^2 + model_fl[*,5]^2))
  mu_fl_i = 1d/((1d - int_fl[*,2])^2 - (int_fl[*,3]^2 + int_fl[*,4]^2))

  plot, mu_fl_m,mu_fl_i,/psym,$
        title='mu FL',xtitle='model value',ytitle='interpreted value'
  oplot,[-500,500],[-500,500]
  device,/close
  set_plot,'X'

  set_plot,'PS'
  !p.font=0
  device,filename='verify_FL_fflex1.eps',$
         xsize=6,ysize=6,/inches,/decomposed,/color,bits_per_pixel=8,/encapsulated
  plot, model_fl[*,6],int_fl[*,5],/psym,$
        title='fflex1 FL',xtitle='model value',ytitle='interpreted value'
  oplot,[-500,500],[-500,500]
  device,/close
  set_plot,'X'

  set_plot,'PS'
  !p.font=0
  device,filename='verify_FL_fflex1.eps',$
         xsize=6,ysize=6,/inches,/decomposed,/color,bits_per_pixel=8,/encapsulated
  plot, model_fl[*,6],int_fl[*,5],/psym,$
        title='fflex1 FL',xtitle='model value',ytitle='interpreted value'
  oplot,[-500,500],[-500,500]
  device,/close
  set_plot,'X'

  set_plot,'PS'
  !p.font=0
  device,filename='verify_FL_fflex2.eps',$
         xsize=6,ysize=6,/inches,/decomposed,/color,bits_per_pixel=8,/encapsulated
  plot, model_fl[*,7],int_fl[*,6],/psym,$
        title='fflex2 FL',xtitle='model value',ytitle='interpreted value'
  oplot,[-500,500],[-500,500]
  device,/close
  set_plot,'X'

  set_plot,'PS'
  !p.font=0
  device,filename='verify_FL_gflex1.eps',$
         xsize=6,ysize=6,/inches,/decomposed,/color,bits_per_pixel=8,/encapsulated
  plot, model_fl[*,8],int_fl[*,7],/psym,$
        title='gflex1 FL',xtitle='model value',ytitle='interpreted value'
  oplot,[-500,500],[-500,500]
  device,/close
  set_plot,'X'

  set_plot,'PS'
  !p.font=0
  device,filename='verify_FL_gflex2.eps',$
         xsize=6,ysize=6,/inches,/decomposed,/color,bits_per_pixel=8,/encapsulated
  plot, model_fl[*,9],int_fl[*,8],/psym,$
        title='gflex2 FL',xtitle='model value',ytitle='interpreted value'
  oplot,[-500,500],[-500,500]
  device,/close
  set_plot,'X'

endif




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


end

