function sis_lensmodel, coords, theta_e=theta_e, polar=polar, kappa=kappa, $
                        ctr=ctr, inc_pos=inc_pos, inc_kappa=inc_kappa,$
                        psipars=psipars,flexionscale=flexionscale

; This function returns the six lensing parameters of reduced shear
; and flexion assuming a SIS model for a lens for each of the
; positions input in COORDS.  COORDS is assumed to be an N by 2 array
; of positions.  It is assumed that only the lens parameters will be
; returned. However, if the COMBO keyword is set, then an N by 8 array
; will be returned with the first two columns containing the positions
; (model[*,0:1]=coords).  If the INC_KAPPA keyword is set, the last
; column will be filled with the convergence values exported into kappa.
;
; The fields returned are those in Schneider and Er 2008 formalism,
; namely g, G1 and G3, all as complex fields.

; Some protection vs bad calls
if not keyword_set(theta_e) then theta_e=0d
if not keyword_set(ctr) then ctr=[0d,0d]
if not keyword_set(flexionscale) then flexionscale=1d3

theta_e=double(theta_e)

dims=size(coords,/dimensions)
if n_elements(dims) eq 1 then coords=transpose(coords)


nobj=(size(coords,/dimensions))[0]

; Set up the outputs
kappa=dblarr(nobj)
model=dblarr(nobj,6)

; Don't enter the loop if nothing is there


if theta_e gt 0d then for i=0,nobj-1 do begin
; Now to the lens model.
   if not keyword_set(polar) then begin
      r=sqrt(total((coords[i,0:1]-ctr)^2,/double))/theta_e
      phi=atan(coords[i,1]-ctr[1],coords[i,0]-ctr[0])
   endif else begin
      r=coords[i,0]/theta_e
      phi=coords[i,1]
   endelse

   kappa[i]=(0.5d)/r

   g=-(1d)/((2d)*r-1d)*[cos((2d)*phi),sin((2d)*phi)]
   
   G1=flexionscale*(2d)*([cos(phi),sin(phi)]/theta_e)*$
      ((1d)/r-(1d))/((2d)*r-1d)^2
   
   G3=flexionscale*(2d)*([cos((3d)*phi),sin((3d)*phi)]/theta_e)*$
      (3d - (1d)/r)/((2d)*r-1d)^2

   model[i,*]=[g,G1,G3]
   
   if keyword_set(psipars) then $
      model[i,*]=(convert_flexpars([dblarr(6),g,G1,G3],/gn_to_psin))[6:11]

endfor

if keyword_set(inc_pos) then begin
   pmodel=dblarr(nobj,8)
   pmodel[*,0:1]=coords
   pmodel[*,2:7]=model
   model=pmodel
endif

if keyword_set(inc_kappa) then begin
   dims=size(model,/dimensions)+[0,1]
   kmodel=dblarr(dims)
   kmodel[*,0:dims[1]-2]=model
   kmodel[*,dims[1]-1]=kappa
   model=kmodel
endif

if nobj eq 1 then model=transpose(model)

return,model

end
