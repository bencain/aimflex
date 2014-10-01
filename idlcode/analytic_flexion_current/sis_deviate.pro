function sis_deviate, scaled_pars, flexions=flexions, flexionerrs=flexionerrs,$
                      pmids=pmids, pranges=pranges


pars=scaled_pars*pranges+pmids

theta_e=pars[0]
rimage=theta_e*pars[1]
phi=pars[2]

model_psi1=-0.25d*theta_e/(2d*rimage^2-rimage*theta_e)*[cos(phi),sin(phi)]
model_psi3= 0.75d*theta_e/(2d*rimage^2-rimage*theta_e)*[cos(3d*phi),sin(3d*phi)]

if not keyword_set(flexionerrs) then flexionerrs=dblarr(4)+1d

dev=(flexions-[model_psi1,model_psi3])/flexionerrs


return,dev

end
