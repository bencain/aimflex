function unlensed_beta, theta, z=z, lens=lens

; This function returns the source plane position given a set lensing
; potential. The pars variable will hold the Z_weight, then the lens
; model parameters. The theta variable is assumed to be a single object
; position. 

  alpha=my_lens(theta[0:1],lens,Z_weight=z,/alpha_only)

  return, theta[0:1]-alpha

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function defl_deviate, theta, beta=beta, z=z, lens=lens, err=err

; This function is the deviate needed for MPFIT to minimize the
; difference between the source position beta and the source position
; implied by the lens and the image position. We use this to find the
; value of theta which matches the input beta.

  dev=(beta - unlensed_beta(theta,z=z,lens=lens))/abs(err[0])

  return, dev

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function find_image_pos, beta, theta_start, err, Z_src, lens, chisq=chisq

  b=[beta[0],beta[1]]
  s=[theta_start[0],theta_start[1]]
  z=Z_src[0]
  e=abs(err[0])

; Set up the minimization stuff
  fname='defl_deviate'

  parinf=replicate({value:0d, fixed:0, limited:[0,0], limits:[-1d3,1d3],$
                    parname:'',mpmaxstep:0d},2)


  fargs={beta:b,z:z,lens:lens,err:e}

  parinf.value=s

  theta=mpfit(fname,functargs=fargs,parinfo=parinf,/quiet,bestnorm=chisq,maxiter=500L)

  return, theta
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function find_mult_images, beta, err, Z, lens, $
                           ntries=ntries, seed=seed,span=span, $
                           theta_starts=theta_starts, nimages=nimages,$
                           chisq=chisq, show=show

; This function tests for multiple images by solving for image
; solutions by solving the lens equation for a large number of input
; starting values.

  if not keyword_set(ntries) then ntries=625L
  if ntries lt 16L then ntries=16L
  if not keyword_set(span) then span=2d
  
  theta_starts=dblarr(ntries,2)
  theta_mins=dblarr(ntries,2)
  chisq_mins=dblarr(ntries)

  theta_starts[*,0]=span*(randomu(seed,ntries,/double) - 0.5d) + beta[0]
  theta_starts[*,1]=span*(randomu(seed,ntries,/double) - 0.5d) + beta[1]

  start=dblarr(2)
  
  for i=0,ntries-1 do begin

     start[0:1]=theta_starts[i,0:1]
     
     theta=find_image_pos(beta,start,err,Z,lens,chisq=X2)
     theta_mins[i,*]=theta[*]
     chisq_mins[i]=X2

  endfor

; Check to see if any of the solutions are actually ok

  conv=where(sqrt(chisq_mins) lt err,nconv) ; These are image positions where beta is reproduced
  
  if nconv gt 0 then begin
     theta_mins=theta_mins[conv,*] 
     chisq_mins=chisq_mins[conv]
  endif else begin
     nimages=0
     chisq=999d
     return,-1
  endelse

; Eliminate doubles
  theta_outs=dblarr(nconv,2)
  chisq_outs=dblarr(nconv)

  nimages=0
  for i=0,nconv-1 do begin
     use=1
     for j=i+1,nconv-1 do begin
        if sqrt(total((theta_mins[i,*]-theta_mins[j,*])^2)) lt err then use=0
     endfor
     if use then begin
        theta_outs[nimages,*]=theta_mins[i,*]
        chisq_outs[nimages] = chisq_mins[i]
        nimages++
     endif
  endfor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  if keyword_set(show) then print,'THIS IS ALL WRONG - show does not work'
     
;     plot,[-span[0]/2d,span[0]/2d,span[0]/2d,-span[0]/2d,-span[0]/2d],$
;          [-span[0]/2d,-span[0]/2d,span[0]/2d,span[0]/2d,-span[0]/2d],$
;          psym=3,/isotropic

;     ang=2d*!dpi*dindgen(100)/99d

; plot the scale of the first lens
;     theta_e = lens[2]
;     oplot,theta_e*cos(ang),theta_e*sin(ang)


; Put in the source and image positions
;     oplot,[beta[0]],[beta[1]],psym=1,color=fsc_color('green'),thick=2,symsize=2
;     oplot,[theta_outs[u,0]],[theta_outs[u,1]],psym=4,color=fsc_color('red'),thick=2,symsize=2

;     print,'nimages     = ',nimages
;     print,'Source ==>> green +'
;     print,'beta[0]     = ',beta[0]
;     print,'beta[1]     = ',beta[1]
;     print,'Image(s) ==>> red diamond(s)'
;     print,'theta[*,0]  = ',theta_outs[u,0]
;     print,'theta[*,1]  = ',theta_outs[u,1]
;     print,'---'    
;  endif
 
  chisq= chisq_outs[0:nimages-1]
  return,theta_outs[0:nimages-1,*]

end
