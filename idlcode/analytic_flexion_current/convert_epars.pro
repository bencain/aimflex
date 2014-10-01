function convert_epars, pars, ellipse_pars=ellipse_pars,$
                        ae_to_pol=ae_to_pol, pol_to_mi=pol_to_mi,$
                        mi_to_ae=mi_to_ae, pol_to_ae=pol_to_ae,$
                        mi_to_pol=mi_to_pol, ae_to_mi=ae_to_mi

; This function converts between parametrizations of ellipse
; parameters. Both consider an elliptical Gaussian profile given by 
; 
; I(x,y)=S0*exp(-0.5 * r^2) where
;
; r^2 = ((x-X_c)*cos(xi)+(y-Y_c)*sin(xi))^2/A^2 + 
;       ((y-Y_c)*cos(xi)-(x-X_c)*sin(xi))^2/(eps*A)^2
;
; Non-ellipse parameters are unchanged by this function.  If the
; ellipse pars are not specified, we assume they are in pars[3:5] and
; are defined by one of the following parametrizations:
;
; AE:  Semimajor/minor axis geometric mean (AB)^0.5, 
;      axis ratio EPS and position angle XI
; POL: alpha=(AB)^0.5, "plus-mode" ellipticity EPLUS and
;      "cross-mode" ellipticity ECROSS (evoking polarization)
;          r^2=(1/alpha^2) * (1/sqrt(1-eplus^2-ecross^2) *
;              ((1-eplus)*(x-X_c)^2 + (1+eplus)*(y-Y_c) -
;               2*ecross*(x-X_c)*(y-Y_c))
; MI:  Inverse x^2 coefficient M1, inverse y^2 coefficient M2 and xy
;      coefficient M3 from 
;          r^2 = ((x-X_c)^2/M1^2 + (y-Y_c)^2/M2 + 2*M3*(x-X_c)*(y-Y_c))
;
; For Sersic profile use, the ellipses are defined the same way, the
; profile is just defined with index n as:
;              I ~ exp(-0.5 * r^(1/n))
;

if not keyword_set(ellipse_pars) then ellipse_pars=[3,4,5]
outpars=pars


if keyword_set(ae_to_pol) then begin
   A=pars[ellipse_pars[0]]
   eps=pars[ellipse_pars[1]]
   xi=pars[ellipse_pars[2]]

   alpha=A*sqrt(eps)

; This is 
;  E = ( Q11 - Q22 + 2iQ12 )/( Q11 + Q22 )
;   eplus= ((1d)-eps^2)*cos((2d)*xi)/((1d)+eps^2)
;   ecross=((1d)-eps^2)*sin((2d)*xi)/((1d)+eps^2)

; This is 
;  E = ( Q11 - Q22 + 2iQ12 )/( Q11 + Q22 +2sqrt( Q11*Q22 - Q12^2 ) )
   eplus= ((1d)-eps)*cos((2d)*xi)/((1d)+eps)
   ecross=((1d)-eps)*sin((2d)*xi)/((1d)+eps)

   outpars[ellipse_pars]=[alpha,eplus,ecross]

endif

if keyword_set(pol_to_mi) then begin
   alpha=pars[ellipse_pars[0]]
   eplus=pars[ellipse_pars[1]]
   ecross=pars[ellipse_pars[2]]
   emag=sqrt(eplus^2+ecross^2)

; This is 
;  E = ( Q11 - Q22 + 2iQ12 )/( Q11 + Q22 )
;   m1=alpha^2*sqrt(1d -emag^2)/(1d -eplus)
;   m2=alpha^2*sqrt(1d -emag^2)/(1d +eplus)
;   m3=-ecross/(alpha^2*sqrt(1d -emag^2))

; This is 
;  E = ( Q11 - Q22 + 2iQ12 )/( Q11 + Q22 +2sqrt( Q11*Q22 - Q12^2 ) )
   m1=alpha^2*(1d - emag^2)/(1d + emag^2 - 2d*eplus)
   m2=alpha^2*(1d - emag^2)/(1d + emag^2 + 2d*eplus)
   m3=-2d*ecross/(alpha^2*(1d -emag^2))

   outpars[ellipse_pars]=[m1,m2,m3]
endif


if keyword_set(pol_to_ae) then begin
   alpha=pars[ellipse_pars[0]]
   eplus=pars[ellipse_pars[1]]
   ecross=pars[ellipse_pars[2]]
   emag=sqrt(eplus^2+ecross^2)

; This is 
;  E = ( Q11 - Q22 + 2iQ12 )/( Q11 + Q22 )
;   eps=sqrt((1d -emag)/(1d +emag))
;   A=alpha/sqrt(eps)
;   xi=0.5d*atan(ecross,eplus)

; This is 
;  E = ( Q11 - Q22 + 2iQ12 )/( Q11 + Q22 +2sqrt( Q11*Q22 - Q12^2 ) )
   eps=(1d - emag)/(1d + emag)
   A=alpha/sqrt(eps)
   xi=0.5d*atan(ecross,eplus)

   outpars[ellipse_pars]=[A,eps,xi]

endif

if keyword_set(mi_to_pol) then begin
   m1=pars[ellipse_pars[0]]
   m2=pars[ellipse_pars[1]]
   m3=pars[ellipse_pars[2]]

; This is 
;  E = ( Q11 - Q22 + 2iQ12 )/( Q11 + Q22 )
;   alpha=(m1*m2/(1d - m1*m2*m3^2))^(0.25d)
;   eplus=(m1-m2)/(m1+m2)
;   ecross=-2d*m1*m2*m3/(m1+m2)

; This is 
;  E = ( Q11 - Q22 + 2iQ12 )/( Q11 + Q22 +2sqrt( Q11*Q22 - Q12^2 ) )
   alpha=(m1*m2/(1d - m1*m2*m3^2))^(0.25d)
   eplus=(m1 - m2)/(m1 + m2 + 2d*sqrt(m1*m2 - (m1*m2*m3)^2)) 
   ecross=-2d*m1*m2*m3/(m1 + m2 + 2d*sqrt(m1*m2 - (m1*m2*m3)^2))

   outpars[ellipse_pars]=[alpha,eplus,ecross]

endif

if keyword_set(ae_to_mi) then $
   outpars=convert_epars(convert_epars(pars,/ae_to_pol),/pol_to_mi)

if keyword_set(mi_to_ae) then $
   outpars=convert_epars(convert_epars(pars,/mi_to_pol),/pol_to_ae)

return, outpars

end
