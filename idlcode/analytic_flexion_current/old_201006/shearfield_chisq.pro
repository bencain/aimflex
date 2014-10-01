function shearfield_chisq, NEW_SHEARS, LENSMODEL_NOW=LENSMODEL_NOW, $
                           WEIGHTS=WEIGHTS, FLEXIONSCALE=FLEXIONSCALE,$
                           CHISQ_DIFF=CHISQ_DIFF, CHISQ_NORM=CHISQ_NORM

; For a set of N galaxies, NEW_SHEARS will be a 2N element array of
; the new shears.  Elements 0 => N-1 are the real components and
; elements N => 2N-1 are the imaginary components.  LENSMODEL_NOW is
; an Nx8 array with the current lens parameters.  Columns 0-1 are the
; x and y positions.  Columns 2-3 are the current real and imaginary
; shear components, respectively.  Columns 4-5 are the current
; 1-flexion and columns 6-7 are the current 3-flexion.  WEIGHTS is an
; NxN array of the pair weightings.

n=(size(lensmodel_now,/dimensions))[0]
if not keyword_set(flexionscale) then flexionscale=1d3

ks=lindgen(2*n^2)
is=(ks mod n^2) mod n
js=long(ks mod n^2)/long(n)
ms=long(ks)/long(n^2)

kpair=ks[where(is ne js,npair)]
ipair=is[kpair]
jpair=js[kpair]
mpair=ms[kpair]

delta_x=dblarr(npair)
delta_y=dblarr(npair)

grad_g_x=dblarr(npair)
grad_g_y=dblarr(npair)

flex_shear_diff=dblarr(npair)
new_shear_diff=dblarr(npair)

ws=dblarr(npair)

delta_x[*]=lensmodel_now[ipair,0]-lensmodel_now[jpair,0]
delta_y[*]=lensmodel_now[ipair,1]-lensmodel_now[jpair,1]

grad_g_x[*]=0.5d*(lensmodel_now[jpair,mpair+6L]+$
                  lensmodel_now[jpair,mpair+4L])/flexionscale
grad_g_y[*]=$
   0.5d*(-1d)^mpair*(lensmodel_now[jpair,((mpair+1) mod 2)+6L]-$
                     lensmodel_now[jpair,((mpair+1) mod 2)+4L])/flexionscale

flex_shear_diff[*]=delta_x*grad_g_x+delta_y*grad_g_y
new_shear_diff[*]=new_shears[ipair+mpair*n]-new_shears[jpair+mpair*n]

ws[*]=weights[ipair,jpair]

chisq_diff=total(ws*(flex_shear_diff-new_shear_diff)^2)

; Add in the part which keeps the shears near zero.
sig_g=0.3d
reg_scale=0d ; This could be changed


chisq_norm=total(reg_scale*(mean(new_shears)/sig_g)^2)


return,chisq_diff+chisq_norm

end
