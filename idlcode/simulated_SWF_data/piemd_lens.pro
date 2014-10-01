function piemd_lens, theta, pars, Z_weight=Z_weight, alpha_only=alpha_only

; This function returns the lensing fields for a non-singular
; isothermal ellipse truncated at some radius.  The output is:
;
;    psi, alpha1, alpha2, kappa, gamma1, gamma2, 
;    fflex1, fflex2, gflex1, gflex2
;
; The input is the position (theta) and the parameters which describe
; the lensing system (pars).  If theta is a 1-D input, then theta[0:1]
; gives the position of the object. If theta is multidimensional, then
; it is taken to be nobj x 2, where nobj is the number of image
; positions. 

; pars gives the details of the lens, in order: center x, y, the Einstein
; radius b for a source at infinity, the softening radius s, the
; truncation radius t, the axis ratio q, and the position angle a.

; The Z_weight keyword gives the scaling from a source at infinity to
; a finite source distance.

  nie_pos_pars=pars[[0,1,2,3,5,6]]
  nie_neg_pars=pars[[0,1,2,4,5,6]]

  dim_t=size(theta,/n_dimensions)

  if dim_t eq 1 then begin
     nobj=1
  endif else begin
     nobj=n_elements(theta[*,0])
  endelse

  if not keyword_set(Z_weight) then begin
     Z=1d 
  endif else begin
     if n_elements(Z_weight) lt nobj then begin
        Z=Z_weight[0] 
     endif else begin
        Z=Z_weight[0:nobj-1]
     endelse
  endelse

  pos_fields=nie_lens(theta,nie_pos_pars,Z_weight=Z,alpha_only=keyword_set(alpha_only))
  neg_fields=nie_lens(theta,nie_neg_pars,Z_weight=Z,alpha_only=keyword_set(alpha_only))

  return,pos_fields-neg_fields

end
