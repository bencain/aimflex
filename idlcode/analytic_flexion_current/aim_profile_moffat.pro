function aim_profile_moffat, mpars, beta

; This function takes in a set of model parameters (MPARS) and a
; complex coordinate array (BETA) and returns the surface brightness
; associated with the appropriate elliptical Moffat profile.  I assume
; that  BETA has been calculated from a center-subtracted image plane
; position (therefore mpars[1:2] are ignored)

I0=(10d)^mpars[0]/(2d*!dpi*mpars[3]^2)
E=dcomplex(mpars[4],mpars[5])

rsq=real_part((1d + E*conj(E))*beta*conj(beta) - beta^2*conj(E) - conj(beta)^2*E)

return,I0*(1d + rsq/mpars[3]^2)^(-mpars[6])

end
