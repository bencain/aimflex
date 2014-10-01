pro paper_data_multiple

  ra=104.5665559d               ; DD
  dec=-55.9434082d              ; DD
  hdrfile='~/idl/pro/simulated_SWF_data/1e0657R_avg_wcs_hdr.fits'

  density=80d                   ; arcmin^-2
  field=3.5d                    ; arcmin


  nreal=10
  nsmall=7
  nbig=3
  names='simulation'+string(lindgen(nreal)+1,format='(I02)')

  q=[replicate(1d,nsmall-1),0.75d,replicate(1d,nbig-1),0.75] ; Axis ratios
  position_angle=randomu(seed,nreal,/double)*!dpi       ; Position angles (in radians)

  fieldsize=[replicate(1d,nsmall),replicate(2d,nbig)]*field

  sigv1=[0.5d, randomu(seed,nreal-1)*0.3d + 0.2d ]
  sigv2=[0.5d, randomu(seed,nreal-1)*0.3d + 0.2d ]


; Set up the central halo
  ctr=[0d,0d]                   ; offset the NIE (in arcmin)
  sigma=1.5d                    ; Velocity Dispersion, in 10^3 km/s
  kappa0=0d                     ; Mass sheet offset
  theta_c=1d-2                  ; Core radius, in arcmin

  epsilon=(1d - q)/(1d + q)     ; Ellipticity

  z_lens=0.3d                   ; Lens redshift
  DLI_DI = ( redshift_to_angdd(1d4) - (1d + z_lens)*redshift_to_angdd(z_lens)/(1d + 1d4) )/redshift_to_angdd(1d4)
                                ; Lensing weight


  te0=0.5d*(sigma)^2  * DLI_DI
  te1=0.5d*(sigv1)^2  * DLI_DI
  te2=0.5d*(sigv2)^2  * DLI_DI

;  print,te0
;  print,te1
;  print,te2

; Make a set of realizations
  for i=nreal-1,nreal-1 do begin

;Substructures
     subok=0

     mainhalo=[ctr[0],ctr[1], $
               sigma,kappa0,theta_c, $
               epsilon[i],position_angle[i],DLI_DI] ;horizontal ellipse

     print,epsilon[i],position_angle[i]

     while not subok do begin

        r=fieldsize[i]*(0.3d*randomu(seed,2) + 0.1d) ; between a maximum radius of 0.4 and a minimum of 0.1 times the field size.
        ; This is 20 - 80% of the way from the center to the edge
        a=2d*!dpi*randomu(seed,2)
        sub1=[r[0]*cos(a[0]),r[0]*sin(a[0]), sigv1[i], 0d, 5d-3, 0d,0d, DLI_DI]
        sub2=[r[1]*cos(a[1]),r[1]*sin(a[1]), sigv2[i], 0d, 5d-3, 0d,0d, DLI_DI]

        pts=2d*!dpi*dindgen(100)/99d
        xhat=cos(pts)
        yhat=sin(pts)
        
        plot,[-0.5,0.5,0.5,-0.5,-0.5]*fieldsize[i],[-0.5,-0.5,0.5,0.5,-0.5]*fieldsize[i],/isotropic
        oplot,[ctr[0]],[ctr[1]],psym=4,symsize=3,color=fsc_color('red'),thick=3
        oplot,[sub1[0]],[sub1[1]],psym=4,symsize=3,color=fsc_color('blue'),thick=3
        oplot,[sub2[0]],[sub2[1]],psym=4,symsize=3,color=fsc_color('blue'),thick=3

        oplot,te0*xhat+ctr[0], te0*yhat+ctr[1]
        oplot,te0*q[i]*xhat+ctr[0], te0*q[i]*yhat+ctr[1]
        oplot,[ctr[0],ctr[0]+te0*cos(position_angle[i])],[ctr[0],ctr[0]+te0*sin(position_angle[i])]

        oplot,te1[i]*xhat+sub1[0], te1[i]*yhat+sub1[1]
        oplot,te2[i]*xhat+sub2[0], te2[i]*yhat+sub2[1]

        read,'Look good? (1=y,0=n)',subok
     endwhile

     nr=1
     nsl=5
  
  
     simdata,names[i],density,fieldsize[i],nr,z_lens=z_lens,seed=seed, n_slsys=nsl, bcg_ra=ra,bcg_dec=dec,$
             lens0=mainhalo,lens1=sub1, lens2=sub2

  endfor

end
