;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION NFW_XI,X, TOL=TOL
; Helper function for NFW lensing from Leonard et al. 2010
  if not keyword_set(tol) then tol=1d-6

  eq1=where(abs(x-1d) lt tol/2d,neq1)

  gt1=where(x ge 1d + tol,ngt1)
  lt1=where(x le 1d - tol,nlt1)

  xiofx=x*0d

  if nlt1 ne 0 then xiofx[lt1]=$
     (2d/sqrt(1d - x[lt1]^2))*atanh(sqrt((1d - x[lt1])/(1d + x[lt1])))

  if ngt1 ne 0 then xiofx[gt1]=$
     (2d/sqrt( x[gt1]^2 - 1d))*atan(sqrt((x[gt1] - 1d)/(x[gt1] - 1d)))

  if neq1 gt 0 then begin
     xiofx[eq1]=0.5d*(nfw_xi(1d - 2d*tol) + nfw_xi(1d + 2d*tol))
  endif

  return,xiofx
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
