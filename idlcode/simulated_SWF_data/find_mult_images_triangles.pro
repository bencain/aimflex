function find_images, beta, Z, lens, center, span, ngrid, $
                      src_dist=src_dist, nimages=nimages, show=show

; This function searches within the limits given to give a new set of
; limits on the image positions based on mapping the triangles back to
; the source plane and overlapping them with the source position.

  xg=span*(dindgen(ngrid^2) mod ngrid)/(ngrid-1d) - span/2d + center[0]
  yg=span*double(lindgen(ngrid^2)/ngrid)/(ngrid-1d) - span/2d + center[1]
  tg=dblarr(ngrid^2,2)
  tg[*,0]=xg
  tg[*,1]=yg

  bg=tg-my_lens(tg,lens,Z_weight=Z,/alpha_only) ; Where do the grid points map

  triangulate,xg,yg,tri         ; Get image plane triangle nodes
  ntri=n_elements(tri[0,*])

  if keyword_set(show) then begin
     ; plot the triangle grid
     minx=min(tg[*,0])
     maxx=max(tg[*,0])
     miny=min(tg[*,1])
     maxy=max(tg[*,1])
     plot,[minx,maxx,maxx,minx,minx],$
          [miny,miny,maxy,maxy,miny],$
          /isotropic
     oplot,[beta[0]],[beta[1]],color=fsc_color('red'),symsize=2
  endif

  tri_bx=tri*0d
  tri_by=tri*0d

  tri_bx[0,*]=bg[tri[0,*],0]    ; Get source plane coords for triangle nodes
  tri_bx[1,*]=bg[tri[1,*],0]
  tri_bx[2,*]=bg[tri[2,*],0]

  tri_by[0,*]=bg[tri[0,*],1]
  tri_by[1,*]=bg[tri[1,*],1]
  tri_by[2,*]=bg[tri[2,*],1]

  dx= -tri_bx + beta[0] ; Get components of vectors from node to source pos
  dy= -tri_by + beta[1]

; Do some cross products.  Note that there are some sign
; problems in Bartelmann's paper.
  d0xd1 = dx[0,*]*dy[1,*] - dy[0,*]*dx[1,*]
  d1xd2 = dx[1,*]*dy[2,*] - dy[1,*]*dx[2,*]
  d2xd0 = dx[2,*]*dy[0,*] - dy[2,*]*dx[0,*]

; Which triangles match up to the source position?
  triok=where(((d0xd1 ge 0d) and (d1xd2 ge 0d) and (d2xd0 ge 0d)) or $
              ((d0xd1 lt 0d) and (d1xd2 lt 0d) and (d2xd0 le 0d)),nimages)

  if keyword_set(show) then begin
     for i=0,nimages-1 do begin
        oplot,xg[[tri[0,triok[i]],tri[1,triok[i]],$
                  tri[2,triok[i]],tri[0,triok[i]]]],$
              yg[[tri[0,triok[i]],tri[1,triok[i]],$
                  tri[2,triok[i]],tri[0,triok[i]]]],$
              color=fsc_color('green')
        oplot,[tri_bx[0,triok[i]],tri_bx[1,triok[i]],$
               tri_bx[1,triok[i]],tri_bx[0,triok[i]]],$
              [tri_by[0,triok[i]],tri_by[1,triok[i]],$
               tri_by[1,triok[i]],tri_by[0,triok[i]]],$
              color=fsc_color('blue')
     endfor

     read,blar
  endif

  if nimages lt 1 then begin
     src_dist=9d99
     return,-1
  endif

; We'll return the image plane position of the centroid of the
; nodes and the distance of the unlensed position of the centroid to
; the source in the source plane.
  src_dist=dblarr(nimages)
  out_pos= dblarr(nimages,2)

  for i=0,nimages-1 do begin
     out_pos[i,0] =mean(xg[tri[*,triok[i]]])
     out_pos[i,1] =mean(yg[tri[*,triok[i]]])
  endfor
  diff=out_pos - my_lens(out_pos,lens,Z_weight=Z,/alpha_only)
  diff[*,0]-=beta[0]
  diff[*,1]-=beta[1]
  src_dist=sqrt(total(diff^2,2))

  return,out_pos

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function find_mult_images_triangles, beta, err, Z, lens, $
                                     span=span, ngrid=ngrid, resolution=resolution,$
                                     nimages=nimages, chisq=chisq,show=show

; This function finds multiple images by breaking the lens plane into
; gridded triangles, then determining which triangles contain the
; source position when mapped into the source plane.  See the paper by
; Bartelmann (2003) in Proceedings Contribution, Gravitational Lensing
; Winter School, Aussois 2003 

  if not keyword_set(ngrid) then ngrid=25L
  if ngrid lt 10L then ngrid=10L
  ngrid = long(ngrid)

  if not keyword_set(span) then span=2d
  if span lt 0.1d then span=0.1d
  span=double(span)
  
  xreg_curr = 0d
  yreg_curr = 0d
  nreg_curr = 1
  span_curr = span
  dist_curr = 99d*abs(err)

; We iterate until the positions collected are all within err of the
; source position in the source plane.

  while max(dist_curr) gt abs(err) do begin
     nreg_new=0
; Loop through the regions, find the images there and append them.
     for i=0,nreg_curr-1 do begin
; Check to see whether there is an image here
        tmp_img=find_images(beta, Z, lens, [xreg_curr[i],yreg_curr[i]], span_curr, ngrid,$
                            src_dist=tmp_dist, nimages=ntmp,show=keyword_set(show))
; If there are new images to add, add them        
        if ntmp gt 0 then begin
           if nreg_new eq 0 then begin
              xreg_new=tmp_img[*,0]
              yreg_new=tmp_img[*,1]
              nreg_new=ntmp
              dist_new=tmp_dist
           endif else begin
              xreg_new=[xreg_new,tmp_img[*,0]]
              yreg_new=[yreg_new,tmp_img[*,1]]
              nreg_new+=ntmp
              dist_new=[dist_new,tmp_dist]
           endelse
        endif

     endfor

; Now change the new to current.
     xreg_curr = xreg_new
     yreg_curr = yreg_new
     nreg_curr = nreg_new
     span_curr *= 4d/(ngrid - 1)
     dist_curr = dist_new

  endwhile

; Now we need to get rid of any images that are at the same position:
  xout=xreg_curr[0]             ;Nab the first image
  yout=yreg_curr[0]
  nout=1
  others=where(((xreg_curr - xout[nout-1])^2 + (yreg_curr - yout[nout-1])^2) gt resolution^2, nothers) ; find images away from the first
  while nothers gt 0 do begin
     xreg_curr=xreg_curr[others] ; pare down to just ones away from the current out list
     yreg_curr=yreg_curr[others]
     xout=[xout,xreg_curr[0]]   ; add the first away image
     yout=[yout,yreg_curr[0]] 
     nout++                     ; update the number of images
     others=where(((xreg_curr - xout[nout-1])^2 + (yreg_curr - yout[nout-1])^2) gt resolution^2,nothers) ; find others
  endwhile

  nimages=nout
  images=dblarr(nimages,2)
  images[*,0]=xout
  images[*,1]=yout

  return,images

end
