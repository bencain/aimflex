function find_mult_images_newton,  beta, err, Z, lens, $
                                   span=span, nguess=nguess, niter=niter,$
                                   nimages=nimages,conv_iter=conv_iter,show=show, nhold=nhold

  if not keyword_set(nhold) then nhold=25L
  if nhold lt 10 then nhold=10L
  nhold=long(nhold)

  theta_n=(randomu(seed,nguess,2)-0.5d)*span
  J=dblarr(2,2)
  
  iter=0
  ivn=dblarr(niter)

  if keyword_set(show) then print,'beta = ',beta[0],beta[1]
  while iter lt niter do begin
; Get the fields at the current points
     fields_n=my_lens(theta_n,lens,Z_weight=Z)
     theta_np1=theta_n*0d

     if keyword_set(show) then begin
        plot,theta_n[*,0],theta_n[*,1],/psym,/iso
        oplot,[beta[0]],[beta[1]],/psym,color=fsc_color('blue')
        wait,0.05
     endif
     
     
     for i=0,nguess-1 do begin

        J[0,0]=1d - fields_n[i,3]-fields_n[i,4] ; Lensing Jacobian
        J[1,1]=1d - fields_n[i,3]+fields_n[i,4]
        J[1,0]= -fields_n[i,5]
        J[0,1]= -fields_n[i,5]

        theta_np1[i,*] = theta_n[i,*] + invert(J)##(reform(beta,[1,2]) - theta_n[i,*]-fields_n[i,1:2])

     endfor

     theta_n=theta_np1

;;;;;;;;;;;
     possible=theta_n

     xout=possible[0,0]
     yout=possible[0,1]
     nimages=1
     
     others=where(sqrt((possible[*,0]-xout[nimages-1])^2 + (possible[*,1]-yout[nimages-1])^2) gt err,nothers)
     
     while nothers gt 0 do begin
        possible=possible[others,*]
        xout=[xout,possible[0,0]]
        yout=[yout,possible[0,1]]
        nimages++
        others=where(sqrt((possible[*,0]-xout[nimages-1])^2 + (possible[*,1]-yout[nimages-1])^2) gt err,nothers)
     endwhile 
     ivn[iter]=nimages
;;;;;;;;;;;     
  

; If nothing changes for a while, break out of the loop
;     if (iter gt 3*nhold) then $
;        if (max(abs(ivn[iter-nhold-1:iter-1]-ivn[iter])) eq 0) then break

     iter++

  endwhile

  if keyword_set(show) then begin
     oplot,theta_n[*,0],theta_n[*,1],/psym,color=fsc_color('green')
     wait,0.5
  endif


;  plot,ivn
  conv_iter=min(where(ivn eq min(ivn)))

  possible=theta_n
  xout=possible[0,0]
  yout=possible[0,1]
  nimages=1

  others=where(sqrt((possible[*,0]-xout[nimages-1])^2 + (possible[*,1]-yout[nimages-1])^2) gt err,nothers)

  while nothers gt 0 do begin
     possible=possible[others,*]
     xout=[xout,possible[0,0]]
     yout=[yout,possible[0,1]]
     nimages++
     others=where(sqrt((possible[*,0]-xout[nimages-1])^2 + (possible[*,1]-yout[nimages-1])^2) gt err,nothers)
  endwhile 

  out=dblarr(nimages,2)
  out[*,0]=xout
  out[*,1]=yout

  return,out

end



