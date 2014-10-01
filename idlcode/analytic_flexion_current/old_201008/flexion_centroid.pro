function FLEXION_CENTROID, MODELPARS 
                                ; N x 12 array of fit parameters 

mps=modelpars
dims=size(mps,/dimensions)
if n_elements(dims) lt 2 then begin
   mps=transpose(modelpars)
   dims=size(mps,/dimensions)
endif
if dims[1] lt 12 then begin
   print, 'Error - bad inputs to FLEXION_CENTROID'
   print, 'Must have dimensions of Nx12 or just 12'
   return,0d
endif

centroid=dblarr(dims[0],2)

alpha=modelpars[*,3]
E=dcomplex(modelpars[*,4],modelpars[*,5])
g=dcomplex(modelpars[*,6],modelpars[*,7])
G1=dcomplex(modelpars[*,8],modelpars[*,9])
G3=dcomplex(modelpars[*,10],modelpars[*,11])

Q0=2d*alpha^2*sqrt(1-E*conj(E))
Q2=alpha^2*sqrt((1d +sqrt(E*conj(E)))/(1d -sqrt(E*conj(E))))*E

b1=(2d*g^2*Q0-2d*g*Q2)/(4d*(1d -g*conj(g)))
b2=(8d*g*Q0-5d*Q2-3d*g^2*conj(Q2))/(4d*(1d -g*conj(g)))
b3=(3d*conj(g)*Q2-2d*(3d +g*conj(g))*Q0+5d*g*Q2)/(4d*(1d -g*conj(g)))
b4=((3d*g*conj(g)-1d)*conj(Q2)-2d*conj(g)*Q0)/(4d*(1d -g*conj(g)))

beta=b1*conj(G3)+b2*conj(G1)+b3*G1+b4*G3

centroid[*,0]=real_part(beta)
centroid[*,1]=imaginary(beta)

return,centroid
end
