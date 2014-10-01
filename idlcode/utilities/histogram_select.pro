function histogram_select, arr, name=name, ptitle=ptitle, n_sel=n_sel, $
                           not_sel=not_sel, n_not_sel=n_not_sel,$
                           wtitle=wtitle,$
                           init_range=init_range, range=range, $
                           prompt=prompt, psfile=psfile, toeps=toeps,$
                           onlyupper=onlyupper, onlylower=onlylower

  
; This function plots ARR as a histogram and interactively selects a
; subset of those elements using input from the user.

if not keyword_set(name) then name='Variable 1'
if not keyword_set(wtitle) then wtitle='Select on '+name
if not keyword_set(ptitle) then ptitle=name+' Histogram'
if not keyword_set(prompt) then prompt=$
   'Select the desired elements by '+name

; Mark which ends we're going to check.
do_upper=1
do_lower=1
if keyword_set(onlyupper) then do_lower=0
if keyword_set(onlylower) then do_upper=0

clear
print,' '
print,'********'
print,prompt
print,'********'

window,/free,title=wtitle

a_min=double(min(arr))
a_max=double(max(arr))

width=a_max-a_min
center=(a_max+a_min)/2d

if width eq 0 then width+=0.01d

init_range=center+width*[-0.6d,0.6d]
binsize=(init_range[1]-init_range[0])/100d

yr=[0d,double(n_elements(arr))]

;Check the plot x-range
rng_ok=0
while not rng_ok do begin
   init_range=center+width*[-0.6d,0.6d]
   plot_hist,arr,binsize=binsize,xrange=init_range,peak=peak,$
             title=ptitle,xtitle=name,$
             charsize=2.0,charthick=2.0,thick=2.0,xthick=2.0,ythick=2.0,$
             xmargin=[15,15],ymargin=[8,8]

   clear
   print,'Current center: '+strcompress(string(center),/remove_all)
   print,'Current width:  '+strcompress(string(width),/remove_all)
   rng_q=get_num_resp(['Check the X-axis range.',$
                       '   0 = ok',$
                       '   1 = change center',$
                       '   2 = change width'],$
                      ok_resp=[0,1,2],default=0)

   if rng_q eq 0 then rng_ok=1
   if rng_q eq 1 then begin
      read,center,prompt='Enter the new center: '
      center=double(center)
   endif
   if rng_q eq 2 then begin
      read,width,prompt='Enter the new width: '
      width=abs(double(width))
   endif

   a_min=center-width/2d
   a_max=center+width/2d
endwhile

;Check the plot y-range
rng_ok=0
while not rng_ok do begin
   plot_hist,arr,binsize=binsize,xrange=init_range,peak=peak,$
             title=ptitle,xtitle=name,$
             charsize=2.0,charthick=2.0,thick=2.0,xthick=2.0,ythick=2.0,$
             xmargin=[15,15],ymargin=[8,8]

   clear
   print,'Current y maximum: '+strcompress(string(yr[1]),/remove_all)
   rng_q=get_num_resp(['Check the Y-axis range.',$
                       '   0 = ok',$
                       '   1 = change maximum'],$
                      ok_resp=[0,1],default=0)

   if rng_q eq 0 then rng_ok=1
   if rng_q eq 1 then begin
      read,ym,prompt='Enter the new maximum: '
      yr[1]=abs(double(ym))
   endif

endwhile

; Check the histogram resolution
bins_ok=0
while not bins_ok do begin
   plot_hist,arr,binsize=binsize,xrange=init_range,peak=peak,$
             binsize=binsize,title=ptitle,xtitle=name,$
             charsize=2.0,charthick=2.0,thick=2.0,xthick=2.0,ythick=2.0,$
             xmargin=[15,15],ymargin=[8,8]

   yr[1]=peak*1.2d

   clear
   print,'Current bin size: '+strcompress(string(binsize),/remove_all)
   bins_q=get_num_resp(['Check the bin size.',$
                        '   0 = ok',$
                        '   1 = double bin size',$
                        '   2 = halve bin size',$
                        '   3 = enter custom'],$
                        ok_resp=[0,1,2,3],default=0)

   if bins_q eq 0 then bins_ok=1
   if bins_q eq 1 then binsize*=2d
   if bins_q eq 2 then binsize/=2d
   if bins_q eq 3 then begin
      read,binsize,prompt='Enter the binsize: '
      binsize=double(binsize)
   endif

endwhile




; Start the outlier selection
ok=0
clear
while (not ok) do begin

   if (do_upper and do_lower) then begin
      sel=where((arr le a_max) and (arr ge a_min),n_sel,$
                 complement=not_sel,ncomplement=n_not_sel)
   endif

   if (do_upper and (not do_lower)) then begin
      sel=where(arr le a_max,n_sel,$
                 complement=not_sel, ncomplement=n_not_sel)
   endif

   
   if ((not do_upper) and do_lower) then begin
      sel=where(arr ge a_min,n_sel,$
                 complement=not_sel,ncomplement=n_not_sel)
   endif

   if ((not do_upper) and (not do_lower)) then begin
      sel=-1
      n_sel=0
      not_sel=arr
      n_not_sel=n_elements(arr)
      return,sel
   endif

   plot_hist,arr[sel],xrange=init_range,yrange=yr,$
             binsize=binsize,title=ptitle,xtitle=name,$
             charsize=2.0,charthick=2.0,thick=2.0,xthick=2.0,ythick=2.0,$
             xmargin=[15,15],ymargin=[8,8],ystyle=1,$
             /xhatch,hatch_size=peak/20d

   if n_not_sel gt 0 then plot_hist,arr[not_sel],/overplot,binsize=binsize

   if do_upper then oplot,[a_max,a_max],[0,2d*peak]
   if do_lower then oplot,[a_min,a_min],[0d,2d*peak]

   clear
; Check if the selection is ok
   ok=get_num_resp(['Cross hatched histogram bars are '+$
                    'selected, open bars are not.',$
                    'Is the selection good?',$
                    '  1 = yes;  0 = no'],default=0)

   if (not ok) then begin

      ; Get new limits
      if do_upper then begin
         read,a_max,prompt='Enter new upper bound: '
         oplot,[a_max,a_max],[0,2d*peak],linestyle=3
      endif

      if do_lower then begin
         read,a_min,prompt='Enter new lower bound: '
         oplot,[a_min,a_min],[0d,2d*peak],linestyle=3
      endif

   endif

endwhile

if keyword_set(psfile) then begin

   setup_ps,psfile
   plot_hist,arr[sel],xrange=init_range,yrange=yr,$
             binsize=binsize,xtitle=name,$
             charsize=3.0,charthick=3.0,thick=3.0,xthick=3.0,ythick=3.0,$
             xmargin=[30,30], ymargin=[15,15],ystyle=1,$
             /xhatch,hatch_size=peak/20d
   if n_not_sel gt 0 then plot_hist,arr[not_sel],/overplot,$
                                    binsize=binsize

   setup_x

   print,'Saved to '+psfile
   if keyword_set(toeps) then begin
      ps2eps,psfile,/delete,delok=delok
      if delok then print,psfile+$
                          ' converted to EPS format and PS file deleted.'
   endif

endif

return,sel

end
