pro step_shearfield_test

; Test out the new shearfield fitter.  First some simple data.

n=100

gal_pos=dblarr(n,2)
shear_now=dblarr(n,2)
flex1_now=dblarr(n,2)
flex3_now=dblarr(n,2)


; we'll offset the galaxies by some set grid spacing
fieldsize=4000d
gal_pos[*,0]=(((dindgen(n) mod (n/10))*10d/n) - 0.5d)*fieldsize
gal_pos[*,1]=((double(lindgen(n)/long(n/10))*10d/n) - 0.5d)*fieldsize

;plot,gal_pos[*,0],gal_pos[*,1],psym=1

; Let the shear be zero to start and the shear gradient will be
; constant across the field for each component.

;d1g1=15d
;d2g1=0d

;d1g2=0d
;d2g2=15d

; This means that g1 should be 0 or  +/-0.3 depending on the sign of
; X, and g2 should be 0 or +/-0.3 depending on the sign of y

;flex1_now[*,0]=d1g1+d2g2
;flex1_now[*,1]=d1g2-d2g1

;flex3_now[*,0]=d1g1-d2g2
;flex3_now[*,1]=d1g2+d2g1

lmodel=sis_lensmodel(gal_pos,theta_e=400d,/inc_pos,kappa=kappa)
sub_lmodel=sis_lensmodel(gal_pos,theta_e=0.2d*400d,kappa=subkappa,$
                         ctr=[700d,0d])

for i=0,n-1 do lmodel[i,2:7]=lens_superposition(lmodel[i,2:7],$
                                                sub_lmodel[i,*],$
                                                kappa,subkappa)

good=where((sqrt(total(lmodel[*,4:5]^2,2)) le 15d) or $
          (sqrt(total(lmodel[*,6:7]^2,2)) le 15d),n)

lmodel=lmodel[good,*]

weights=dblarr(n,n)
smoothing=fieldsize*10d/double(n)
for i=0,n-1 do for j=0,n-1 do begin
   dist=sqrt(total((gal_pos[i,*]-gal_pos[j,*])^2))/smoothing
   weights[i,j]=exp(-0.5d*dist^2)
endfor
for i=0,n-1 do weights[i,i]=0d
weights/=total(weights)

;lmodel=dblarr(n,8)
;lmodel[*,0:1]=gal_pos
;lmodel[*,2:3]=shear_now
;lmodel[*,4:5]=flex1_now
;lmodel[*,6:7]=flex3_now

;print,shearfield_deviate(dblarr(2*n),lensmodel_now=lmodel,weights=weights)

;plot_field,lmodel,spin=2
;print,n
;print,size(lmodel,/dimensions)

new_shear=mcmc_shearfield_update(lmodel, weights, errs=shear_errors, $
                                 chisq=chisq)
                                 
new_lmodel=lmodel
new_lmodel[*,2:3]=new_shear

diff_lmodel=new_lmodel
diff_lmodel[*,2:3]-=lmodel[*,2:3]

plot_field,diff_lmodel,spin=2,/autoscale
window,/free
plot_hist,diff_lmodel[*,2:3]
;print,gal_pos
;print,new_shear
;print,shear_errors

end

