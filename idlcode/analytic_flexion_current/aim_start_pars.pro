function aim_start_pars, DATAIMG, BKG, WIN, $
                         NSIGMA=NSIGMA, SIGMA=SIGMA, $
                         SEED=SEED,$
                         ELIMIT=ELIMIT,$
                         SERSIC=SERSIC, PSEUDOGAUSS=PSEUDOGAUSS, $
                         MOFFAT=MOFFAT
  
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

; Find the pixels to use, set the others to zero.
test=dataimg*win-bkg*win

if not keyword_set(sigma) then sigma=stddev(test)
if not keyword_set(nsigma) then nsigma=0d
nsigma=max([nsigma,0d])

low_pix=where(test le nsigma*sigma,nlo,$
              complement=hi_pix,ncomplement=nhi)
if nlo gt 0 then test[low_pix]=0d 
if nhi eq 0d then return,[-5d,0d,0d,2d,dblarr(8)]

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

logN0=alog10(total(test[hi_pix]))
if keyword_set(sersic) then begin
   out=[logN0,Xc,Yc,m1,m2,m3,$
        2d,$ ; Start at a Gaussian
        dblarr(6)] 
endif else if keyword_set(pseudogauss) then begin
   out=[logN0,Xc,Yc,m1,m2,m3,$
        0.9d + 0.2d*randomu(seed,2),$ ; start near a gaussian, but randomized
        dblarr(6)]              
endif else if keyword_set(moffat) then begin
   out=[logN0,Xc,Yc,m1,m2,m3,$
        1.5d,$ ; Start with this, until I get a better sense.
        dblarr(6)]
endif else out=[logN0,Xc,Yc,m1,m2,m3,dblarr(6)]

out=convert_epars(out,/mi_to_pol)

; Rescale the ellipticity if it is over a given limit.
if keyword_set(elimit) then begin
   if elimit le 0d then elimit=1d
   emag=sqrt(total(out[4:5]^2))
   if emag gt elimit then out[4:5]*=0.75d*elimit/emag
endif

; The parameters are:
;        logN0, Xc, Yc, alpha, E+, Ex,
;        [n,] or [k1,k2,]
;        g1, g2, psi11, psi12, psi31, psi32


return,out

end
