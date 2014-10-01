function aim_centroid, g, psi1, psi3, alpha, E

; This function returns the centroid offset from Schneider & Er 2008
; given the apparent ellipse properties and the lensing fields.  g,
; psi1, psi3 and E are assumed to be complex.

Q2=2d*(alpha^2/(1d - abs(E)^2))*(2d*E)
Q0=2d*(alpha^2/(1d - abs(E)^2))*(1d + abs(E)^2)

G1=4d*(psi1 + g*conj(psi1))
G3=4d*(psi3 + g*psi1)

beta=( (3d*conj(g)*G1 - 5d*conj(G1) - 2d*g*conj(G3) )*Q2/$
       (4d*(1d - abs(g)^2)) ) + $
     ( (4d*g*conj(G1) + g^2*conj(G3) - conj(g)*G3 - (3d + abs(g)^2)*G1)*Q0/$
       (2d*(1d - abs(g)^2)) ) + $
     ( (5d*g*G1 - 3d*g^2*conj(G1) - (1d - 3d*abs(g)^2)*G3)*conj(Q2)/$
       (4d*(1d - abs(g)^2)) )

return,beta

end
