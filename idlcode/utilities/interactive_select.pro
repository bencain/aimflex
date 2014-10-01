function interactive_select, arr, brr, aname=aname, bname=bname, $
                             ptitle=ptitle, n_sel=n_sel, not_sel=not_sel,$
                             wtitle=wtitle, n_not_sel=n_not_sel,$
                             ainit_range=ainit_range,$
                             binit_range=binit_range,$
                             a_range=a_range, b_range=b_range, $
                             prompt=prompt, psfile=psfile

  
; This function plots ARR vs BRR and interactively selects a subset of
; those elements using input from the user.

if not keyword_set(aname) then aname='Variable 1'
if not keyword_set(bname) then bname='Variable 2'
if not keyword_set(wtitle) then wtitle='Select in '+$
                                       aname+' vs '+bname+' space.'
if not keyword_set(ptitle) then ptitle=aname+' vs '+bname+'.'
if not keyword_set(prompt) then prompt='Select the desired elements in '+$
                                       aname+' vs '+bname+' space.'

print,' '
print,'********'
print,prompt
print,'********'


amin=double(min(arr))
amax=double(max(arr))

bmin=double(min(brr))
bmax=double(max(brr))

awidth=amax-amin
acenter=(amax+amin)/2d

bwidth=bmax-bmin
bcenter=(bmax+bmin)/2d

if awidth eq 0 then awidth+=0.01d
if bwidth eq 0 then bwidth+=0.01d

ainit_range=acenter+awidth*[-0.6d,0.6d]
binit_range=bcenter+bwidth*[-0.6d,0.6d]

ok=0
aok=0
bok=0

window,xsize=500,ysize=500,/free,title=wtitle

; Start the outlier selection
while (not ok) do begin

; Set the limits for ARR

   while (not aok) do begin

      good=where((arr le amax) and (arr ge amin) and $
                 (brr le bmax) and (brr ge bmin),n_sel,$
                 complement=not_sel,ncomplement=n_not_sel)
      
   ; Recalculate the widths and centers
      awidth=amax-amin
      acenter=(amax+amin)/2d
      
      bwidth=bmax-bmin
      bcenter=(bmax+bmin)/2d

      if awidth eq 0 then awidth+=0.01d
      if bwidth eq 0 then bwidth+=0.01d

      a_range=acenter+(0.55d)*awidth*[-1d,1d]
      b_range=bcenter+(0.55d)*bwidth*[-1d,1d]

      if n_sel gt 0 then begin
         plot,arr[good],brr[good],psym=1,xtitle=aname, ytitle=bname,$
              xrange=ainit_range, yrange=binit_range, title=ptitle
         if n_not_sel gt 0L then oplot,arr[not_sel],brr[not_sel],psym=7
      endif else begin
         plot,arr[not_sel],brr[not_sel],psym=7,xtitle=aname, ytitle=bname,$
              xrange=ainit_range, yrange=binit_range, title=ptitle
      endelse

      oplot,[amin,amax,amax,amin,amin],$
            [bmin,bmin,bmax,bmax,bmin]

   ; Check if the current limits for ARR are ok
      print,'Current center for '+aname+': '+$
            strcompress(string(acenter),/remove_all)
      print,'Current width for '+aname+':  '+$
            strcompress(string(awidth),/remove_all)
      a_q=get_num_resp(['Check current range for '+aname,$
                        '   0 = ok',$
                        '   1 = change center',$
                        '   2 = change width'],$
                       ok_resp=[0d,1d,2d],default=0)

      if a_q eq 0 then aok=1
      if a_q eq 1 then begin
         read,acenter,prompt='Enter the new center: '
         center=double(acenter)
      endif
      if a_q eq 2 then begin
         read,awidth,prompt='Enter the new width: '
         awidth=abs(double(awidth))
      endif
; Reset the limits
      amin=acenter-awidth/2d
      amax=acenter+awidth/2d

   endwhile

; Set the limits for BRR

   while (not bok) do begin

      good=where((arr le amax) and (arr ge amin) and $
                 (brr le bmax) and (brr ge bmin),n_sel,$
                 complement=not_sel,ncomplement=n_not_sel)
      
   ; Recalculate the widths and centers
      awidth=amax-amin
      acenter=(amax+amin)/2d
      
      bwidth=bmax-bmin
      bcenter=(bmax+bmin)/2d

      if awidth eq 0 then awidth+=0.01d
      if bwidth eq 0 then bwidth+=0.01d

      a_range=acenter+0.5d*awidth*[-1d,1d]
      b_range=bcenter+0.5d*bwidth*[-1d,1d]

      if n_sel gt 0 then begin
         plot,arr[good],brr[good],psym=1,xtitle=aname, ytitle=bname,$
              xrange=ainit_range, yrange=binit_range, title=ptitle
         if n_not_sel gt 0L then oplot,arr[not_sel],brr[not_sel],psym=7
      endif else begin
         plot,arr[not_sel],brr[not_sel],psym=7,xtitle=aname, ytitle=bname,$
              xrange=ainit_range, yrange=binit_range, title=ptitle
      endelse
      
      oplot,[amin,amax,amax,amin,amin],$
            [bmin,bmin,bmax,bmax,bmin]

   ; Check if the current limits for BRR are ok
      print,'Current center for '+bname+': '+$
            strcompress(string(bcenter),/remove_all)
      print,'Current width for '+bname+':  '+$
            strcompress(string(bwidth),/remove_all)
      b_q=get_num_resp(['Check current range for '+bname,$
                        '   0 = ok',$
                        '   1 = change center',$
                        '   2 = change width'],$
                       ok_resp=[0,1,2],default=0)

      if b_q eq 0 then bok=1
      if b_q eq 1 then begin
         read,bcenter,prompt='Enter the new center: '
         center=double(bcenter)
      endif
      if b_q eq 2 then begin
         read,bwidth,prompt='Enter the new width: '
         awidth=abs(double(bwidth))
      endif
; Reset the limits
      bmin=bcenter-bwidth/2d
      bmax=bcenter+bwidth/2d

   endwhile

; See if this is ok

   ok_q=get_num_resp(['Check the current ranges',$
                      '   0 = ok',$
                      '   1 = change '+aname,$
                      '   2 = change '+bname,$
                      '   3 = change both'],$
                     ok_resp=[0,1,2,3],default=3)
   
   if ok_q eq 0 then ok=1
   if ok_q eq 1 then begin
      ok=0
      aok=0
      bok=1
   endif
   if ok_q eq 2 then begin
      ok=0
      aok=1
      bok=0
   endif
   if ok_q eq 3 then begin
      ok=0
      aok=0
      bok=0
   endif

endwhile


if keyword_set(psfile) then begin
   setup_ps,psfile

   plot,arr[good],brr[good],psym=1,xtitle=aname, ytitle=bname,$
        xrange=ainit_range, yrange=binit_range,title=ptitle,$
        thick=2,charsize=2,charthick=2
   if n_not_sel gt 0L then oplot,arr[not_sel],brr[not_sel],psym=7
   
   oplot,[amin,amax,amax,amin,amin],$
         [bmin,bmin,bmax,bmax,bmin]
   setup_x
endif

return,good

end
