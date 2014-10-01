function rp_start_pars, DATAIMG, BKG, $
                        NSIGMA=NSIGMA, SIGMA=SIGMA, $
                        SEED=SEED,$
                        ELIMIT=ELIMIT,$
                        SERSIC=SERSIC
  
; This function chooses data-motivated starting values for the fitting
; parameters.  We take only the part of the image which is greater
; than NSIGMA standard deviations of the background level above that
; background level (BKG) assuming Poisson errors in the background.
; Moments are calculated about the center of the image using
; CALC_MOMENTS and inverts those to give approximate starting points.
; If NSIGMA is not set above, then it is set to zero.  We assume that
; the center of DATAIMG is the starting point for the center of the
; unlensed galaxy image.
;
; Note that the starting parameters are returned in the POL format.

if not keyword_set(nsigma) then nsigma=0d
if not keyword_set(sigma) then sigma=stddev(dataimg)

nsigma=max([nsigma,0d])

; Find the pixels to use, set the others to zero.
test=dataimg-bkg
low_pix=where(test lt nsigma*sigma,nlo,complement=hi_pix,ncomplement=nhi)
if nlo gt 0 then test[low_pix]=0d

ratio=(1d)-double(n_elements(low_pix))/double(n_elements(dataimg))

; Make the moments
N=calc_moments(test,0,0)
Xc=calc_moments(test,1,0)/N
Yc=calc_moments(test,0,1)/N
Q11=calc_moments(test,2,0)/N  ; Q11=m1/(1 - m1*m2*m3^2)
Q12=calc_moments(test,1,1)/N  ; Q22=m2/(1 - m1*m2*m3^2)
Q22=calc_moments(test,0,2)/N  ; Q12=-m1*m2*m3/(1 - m1*m2*m3^2)

; Invert the moments to give parameters.
m1=Q11-Q12^2/Q22
m2=Q22-Q12^2/Q11
m3=-Q12/(Q11*Q22-Q12^2)

logN0=alog10(total(dataimg[hi_pix],/double)-double(nhi)*bkg)

if keyword_set(sersic) then out=[logN0,Xc,Yc,m1,m2,m3,0.5d,dblarr(6)] else $
   out=[logN0,Xc,Yc,m1,m2,m3,dblarr(6)]

out=convert_epars(out,/mi_to_pol)

; Rescale the ellipticity if it is over a given limit.
if keyword_set(elimit) then begin
   if elimit le 0d then elimit=1d
   emag=sqrt(total(out[4:5]^2))
   if emag gt elimit then out[4:5]*=0.75d*elimit/emag
endif

; The parameters are:
;        logN0, Xc, Yc, alpha, E+, Ex,
;        [n,]
;        g1, g2, psi11, psi12, psi31, psi32

return,out

end
