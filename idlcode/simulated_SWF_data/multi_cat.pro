pro multi_cat

  ra=104.5665559d               ; DD
  dec=-55.9434082d              ; DD
  hdrfile='~/idl/pro/simulated_SWF_data/1e0657R_avg_wcs_hdr.fits'

  density=80d                  ; arcmin^-2
  field=3.5d                    ; arcmin

; Set up the NIE lens data
  ctr=[0d,0d]                   ; offset the NIE (in arcmin)
  sigma=1.5d                    ; Velocity Dispersion, in 10^3 km/s
  kappa0=0d                     ; Mass sheet offset
  theta_c=1d-2                  ; Core radius, in arcmin
  q=0.75d
  epsilon=(1d - q)/(1d + q)     ; Ellipticity
  position_angle=0d             ; Position angle (in radians)
  z_lens=0.3d                  ; Lens redshift
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
  subok=0
  while not subok do begin

     r=field*0.5d*(randomu(seed,6)*0.5d + 0.25d)
     a=2d*!dpi*randomu(seed,6)
     sub1=[r[0]*cos(a[0]),r[0]*sin(a[0]), 0.5d, 0d, 5d-3, 0d,0d, DLI_DI]
     sub2=[r[1]*cos(a[1]),r[1]*sin(a[1]), 0.5d, 0d, 5d-3, 0d,0d, DLI_DI]
     sub3=[r[2]*cos(a[2]),r[2]*sin(a[2]), 0.5d, 0d, 5d-3, 0d,0d, DLI_DI]
     sub4=[r[3]*cos(a[3]),r[3]*sin(a[3]), 0.5d, 0d, 5d-3, 0d,0d, DLI_DI]
     sub5=[r[4]*cos(a[4]),r[4]*sin(a[4]), 0.5d, 0d, 5d-3, 0d,0d, DLI_DI]
     sub6=[r[5]*cos(a[5]),r[5]*sin(a[5]), 0.5d, 0d, 5d-3, 0d,0d, DLI_DI]

     plot,[-0.5,0.5,0.5,-0.5,-0.5]*field,[-0.5,-0.5,0.5,0.5,-0.5]*field,/isotropic
     oplot,[0d],[0d],psym=4,symsize=3,color=fsc_color('red'),thick=3
     oplot,[sub1[0]],[sub1[1]],psym=4,symsize=3,color=fsc_color('blue'),thick=3
     oplot,[sub2[0]],[sub2[1]],psym=4,symsize=3,color=fsc_color('blue'),thick=3

     read,'Look good? (1=y,0=n)',subok
  endwhile

  nr=1
  nsl=5
  
  
  simdata,'real3',density,field,nr,z_lens=z_lens,seed=seed, n_slsys=nsl, bcg_ra=ra,bcg_dec=dec,$
          lens0=nis;,lens1=sub1, lens2=sub2

end
