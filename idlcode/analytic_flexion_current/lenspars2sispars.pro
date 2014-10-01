FUNCTION LENSPARS2SISPARS, LENSPARS

g=   sqrt(total(lenspars[*,0:1]^2,2))
psi1=sqrt(total(lenspars[*,2:3]^2,2))
psi3=sqrt(total(lenspars[*,4:5]^2,2))

phi= atan(lenspars[*,3],lenspars[*,2])

r= 0.75d*g/psi3

x=(1d + 1d/g)/2d

theta_e=r/x

sispars=dblarr(n_elements(g),3)
sispars[*,0]=r
sispars[*,1]=phi
sispars[*,2]=theta_e


return,sispars

END

