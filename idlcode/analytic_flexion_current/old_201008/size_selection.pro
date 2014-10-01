function size_selection, mag_arr, col_arr, size_arr, $
                         stars=stars, removed=removed


; This function plots color-size and magnitude-size diagrams and takes
; user input to make color/magnitude/size cuts on a galaxy catalog.

if not keyword_set(nsigma) then nsigma=2d

mag_min=min(mag_arr)
mag_max=max(mag_arr)
col_min=min(col_arr)
col_max=max(col_arr)
size_min=min(size_arr)
size_max=max(size_arr)

ok=0
ok1=0
ok2=0
ok3=0

width=mag_max-mag_min
center=(mag_max+mag_min)/2d
init_mrange=[center-(1.1d)*width/2d,center+(1.1d)*width/2d]
width=col_max-col_min
center=(col_max+col_min)/2d
init_crange=[center-(1.1d)*width/2d,center+(1.1d)*width/2d]
width=size_max-size_min
center=(size_max+size_min)/2d
init_srange=[center-(1.1d)*width/2d,center+(1.1d)*width/2d]

window,xsize=500,ysize=800,/free
!p.multi=[0,1,2]

; Loop through until the limits are well set on the magnitude and
; color (to cut out any spurious outliers).
while (not ok) do begin
   
   good=where((mag_arr le mag_max) and (mag_arr ge mag_min) and $
              (col_arr le col_max) and (col_arr ge col_min) and $
              (size_arr le size_max) and (size_arr ge size_min),$
              complement=removed,ncomplement=nbad)

   mwidth=mag_max-mag_min
   mcenter=(mag_max+mag_min)/2d
   magrange=[mcenter-(1.1d)*mwidth/2d,mcenter+(1.1d)*mwidth/2d]

   cwidth=col_max-col_min
   ccenter=(col_max+col_min)/2d
   colrange=[ccenter-(1.1d)*cwidth/2d,ccenter+(1.1d)*cwidth/2d]

   swidth=size_max-size_min
   scenter=(size_max+size_min)/2d
   sizerange=[scenter-(1.1d)*swidth/2d,scenter+(1.1d)*swidth/2d]

; Plot size vs color
   plot,size_arr[good],col_arr[good],psym=1,xtitle='Size (A)',$
        ytitle='Color',xrange=init_srange,yrange=init_crange
   if nbad gt 0 then oplot,size_arr[removed],col_arr[removed],psym=7

   oplot,[size_min,size_max,size_max,size_min,size_min],$
         [col_min,col_min,col_max,col_max,col_min]

; Plot size vs magnitude
   plot,size_arr[good],mag_arr[good],psym=1,xtitle='Size (A)',$
        ytitle='Magnitude',xrange=init_srange,yrange=init_mrange
   if nbad gt 0 then oplot,size_arr[removed],mag_arr[removed],psym=7

   oplot,[size_min,size_max,size_max,size_min,size_min],$
         [mag_min,mag_min,mag_max,mag_max,mag_min]

; Check on the limits for selection.
   if not ok3 then begin
      print,'Current size limits:'
      print,size_min,size_max
      repeat begin
         read,'Are the overall size limits ok (1=yes, 0=no)?  ',ok3 
      endrep until ((ok3 eq 0) or (ok3 eq 1))
   endif
   if not ok3 then begin
      read,'Enter the minimum size for the catalog:  ',size_min
      read,'Enter the maximum size for the catalog:  ',size_max
   endif

   if not ok1 then begin
      print,'Current magnitude limits:'
      print,mag_min,mag_max
      repeat begin
         read,'Are the overall magnitude limits ok (1=yes, 0=no)?  ',ok1 
      endrep until ((ok1 eq 0) or (ok1 eq 1))
   endif
   if not ok1 then begin
      read,'Enter the minimum magnitude for the catalog:  ',mag_min
      read,'Enter the maximum magnitude for the catalog:  ',mag_max
   endif

   if not ok2 then begin
      print,'Current color limits:'
      print,col_min,col_max
      repeat begin
         read,'Are the overall color limits ok (1=yes, 0=no)?  ',ok2 
      endrep until ((ok2 eq 0) or (ok2 eq 1))
   endif
   if not ok2 then begin
      read,'Enter the minimum color for the catalog:  ',col_min
      read,'Enter the maximum color for the catalog:  ',col_max
   endif

   ok=(ok1 and ok2 and ok3)

endwhile

; Now select out the region for red sequence fitting.

ok=0
ok1=0
ok2=0
ok3=0

star_mag_min=mag_min
star_mag_max=mag_max
star_col_min=col_min
star_col_max=col_max
star_size_min=size_min
star_size_max=size_max

while (not ok) do begin
   
   plot,size_arr[good],col_arr[good],psym=2,xrange=sizerange,yrange=colrange,$
        xtitle='Size',ytitle='Color'
   oplot,[star_size_min,star_size_max,$
          star_size_max,star_size_min,star_size_min],$
         [star_col_min,star_col_min,star_col_max,star_col_max,star_col_min]

   plot,size_arr[good],mag_arr[good],psym=2,xrange=sizerange,yrange=magrange,$
        xtitle='Size',ytitle='Magnitude'
   oplot,[star_size_min,star_size_max,$
          star_size_max,star_size_min,star_size_min],$
         [star_mag_min,star_mag_min,star_mag_max,star_mag_max,star_mag_min]
   
   if not ok3 then begin
      print,'Current star size limits:'
      print,star_size_min,star_size_max
      repeat begin 
         read,'Are these limits good (1=yes, 0=no)?  ',ok3 
      endrep until ((ok3 eq 0) or (ok3 eq 1))
   endif
   if not ok3 then begin
      read,'Enter the minimum size for star selection:  ',$
           star_size_min
      read,'Enter the maximum size for star selection:  ',$
           star_size_max
   endif

   if not ok1 then begin
      print,'Current star magnitude limits:'
      print,star_mag_min,star_mag_max
      repeat begin 
         read,'Are these limits good (1=yes, 0=no)?  ',ok1 
      endrep until ((ok1 eq 0) or (ok1 eq 1))
   endif
   if not ok1 then begin
      read,'Enter the minimum magnitude for star selection:  ',$
           star_mag_min
      read,'Enter the maximum magnitude for star selection:  ',$
           star_mag_max
   endif

   if not ok2 then begin
      print,'Current star color limits:'
      print,star_col_min,star_col_max
      repeat begin
         read,'Are these limits good (1=yes, 0=no)?  ',ok2 
      endrep until ((ok2 eq 0) or (ok2 eq 1))
   endif
   if not ok2 then begin
      read,'Enter the minimum color for star selection:  ',star_col_min
      read,'Enter the maximum color for star selection:  ',star_col_max
   endif

   ok=(ok1 and ok2 and ok3)

endwhile

; Plot the two populations
stars=good[where((mag_arr[good] le star_mag_max) and $
                   (mag_arr[good] ge star_mag_min) and $
                   (col_arr[good] le star_col_max) and $
                   (col_arr[good] ge star_col_min) and $
                   (size_arr[good] le star_size_max) and $
                   (size_arr[good] ge star_size_min))]

not_stars=good[where((mag_arr[good] gt star_mag_max) or $
                     (mag_arr[good] lt star_mag_min) or $
                     (col_arr[good] gt star_col_max) or $
                     (col_arr[good] lt star_col_min) or $
                     (size_arr[good] gt star_size_max) or $
                     (size_arr[good] lt star_size_min))]

plot,size_arr[stars],col_arr[stars],psym=2,xrange=sizerange,yrange=colrange,$
     xtitle='Size',ytitle='Color'
oplot,size_arr[not_stars],col_arr[not_stars],psym=4
oplot,[star_size_min,star_size_max,star_size_max,$
       star_size_min,star_size_min],$
      [star_col_min,star_col_min,star_col_max,star_col_max,star_col_min]

plot,size_arr[stars],mag_arr[stars],psym=2,xrange=sizerange,yrange=magrange,$
     xtitle='Size',ytitle='Magnitude'
oplot,size_arr[not_stars],mag_arr[not_stars],psym=4
oplot,[star_size_min,star_size_max,star_size_max,$
       star_size_min,star_size_min],$
      [star_mag_min,star_mag_min,star_mag_max,star_mag_max,star_mag_min]

!p.multi=0

return,not_stars

end




