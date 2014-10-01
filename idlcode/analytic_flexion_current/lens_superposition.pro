function lens_superposition, lenspars_a, lenspars_b, kappa_a, kappa_b, $
                             psipars=psipars

; This function returns the lens parameters of two superimposed
; lenses.  lensparsN=[g1N, g2N, G1_1N, G1_2N, G3_1N, G3_2N].  This
; requires a value for the two convergences.

if keyword_set(psipars) then begin
   lps_a=convert_flexpars(lenspars_a,/psin_to_gn)
   lps_b=convert_flexpars(lenspars_b,/psin_to_gn)
endif else begin
   lps_a=lenspars_a
   lps_b=lenspars_b
endelse

; Get the individual fields
g_a=dcomplex(lps_a[0],lps_a[1])
G1_a=dcomplex(lps_a[2],lps_a[3])
G3_a=dcomplex(lps_a[4],lps_a[5])

g_b=dcomplex(lps_b[0],lps_b[1])
G1_b=dcomplex(lps_b[2],lps_b[3])
G3_b=dcomplex(lps_b[4],lps_b[5])

; Convert back to the physical fields
gamma_a=((1d)-kappa_a)*g_a
fflex_a=(((1d)-kappa_a)/((1d)-g_a*conj(g_a)))*(G1_a-g_a*conj(G1_a))
gflex_a=((1d)-kappa_a)*G3_a-$
        (((1d)-kappa_a)/((1d)-g_a*conj(g_a)))*g_a*(G1_a-g_a*conj(G1_a))

gamma_b=((1d)-kappa_b)*g_b
fflex_b=(((1d)-kappa_b)/((1d)-g_b*conj(g_b)))*(G1_b-g_b*conj(G1_b))
gflex_b=((1d)-kappa_b)*G3_b-$
        (((1d)-kappa_b)/((1d)-g_b*conj(g_b)))*g_b*(G1_b-g_b*conj(G1_b))

; Sum the linear fields
kappa_tot=kappa_a+kappa_b
gamma_tot=gamma_a+gamma_b
fflex_tot=fflex_a+fflex_b
gflex_tot=gflex_a+gflex_b

; Convert back to measurable reduced shear and flexion
g_tot=gamma_tot/((1d)-kappa_tot)
G1_tot=(fflex_tot+g_tot*conj(fflex_tot))/((1d)-kappa_tot)
G3_tot=(gflex_tot+g_tot*fflex_tot)/((1d)-kappa_tot)

lenspars_tot=[real_part(g_tot),imaginary(g_tot),$
              real_part(G1_tot),imaginary(G1_tot),$
              real_part(G3_tot),imaginary(G3_tot)]

if keyword_set(psipars) then $
   lenspars_tot=convert_flexpars(lenspars_tot,/gn_to_psin)

return, lenspars_tot

end
