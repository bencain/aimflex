pro sfm_inp,catfile,sfm,inp,glim=glim,sispars=sispars, $
            fixcentroid=fixcentroid, sersic=sersic, $
            flexlim=flexlim

if not keyword_set(glim) then glim=0.5d
if not keyword_set(flexlim) then flexlim=100d

cat=read_data(catfile)

catsize=size(cat,/dimensions)

if keyword_set(sersic) then npar=13 else npar=12

sfm=dblarr(catsize[0],npar+2)
inp=dblarr(catsize[0],npar)

; 1 = fixed, with a start value
; 2 = with a start value, but free
; Fix the shear, start the flexions
sfm[*,npar-4:npar-3]=1
sfm[*,npar-2:npar+1]=2

;  Get the SIS model lens fields, if parameters given
if not keyword_set(sispars) then begin
   lm=peng_powerlaw_new(cat[*,1:2])
endif else begin
   if n_elements(sispars) ne 3 then lm=peng_powerlaw_new(cat[*,1:2]) else $
      lm=sis_lensmodel(cat[*,1:2],theta_e=sispars[0],ctr=sispars[1:2],$
                       /inc_pos)
endelse

; Put in the starting lens parameters
inp[*,npar-6:npar-1]=lm[*,2:7]

; Fix the center if desired (probably shouldn't do this too much)
if keyword_set(fixcentroid) then sfm[*,3:4]=1

; Set the totals for the started and fixed parameters
for i=0,catsize[0]-1 do begin
   sfm[i,0]=total(sfm[i,2:npar+1] eq 2)
   sfm[i,1]=total(sfm[i,2:npar+1] eq 1)
endfor

; Find the places where the shear input is too large
bad=where(sqrt(total(lm[*,2:3]^2,2)) gt glim,nbad)

; Set bad shear to zero, but fixed.
if nbad gt 0 then inp[bad,npar-6:npar-5]=0d

; find the bad flexions and set them to zero.
bad=where((sqrt(total(lm[*,4:5]^2,2)) gt flexlim) or $
          (sqrt(total(lm[*,6:7]^2,2)) gt flexlim),nbad)
if nbad gt 0 then inp[bad,npar-4:npar-1]=0d

end
