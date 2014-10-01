pro forloop_status,i,n,wnum,label=label, delete=delete

; This procedure plots a status bar for a for loop with some user
; supplied info.

if ((n_elements(wnum) eq 0) or  (!d.window eq -1)) then begin
   window,xsize=450,ysize=150,title='Window opened '+systime(),/free
   wnum=!d.window
endif else begin
   wset,wnum
endelse

if keyword_set(delete) then wdelete,wnum else begin

   if not keyword_set(label) then label='*** FOR-LOOP STATUS ***'

; How far along are we
   pct_done=(100d)*double(i)/double(n)

; Graph points
   pcts=dindgen(101)
   ys=dblarr(101)

   done=where(pcts le pct_done)

   plot,pcts,ys,xrange=[-0.9,100.9],yrange=[-1,1],psym=3,title=label,$
        xtitle='On '+strcompress(string(i+1)+'/'+string(n),/remove_all)+$
        ' ('+strcompress(string(pct_done,format='(D4.1)'),/remove_all)+$
        '% Done)',ystyle=4,xstyle=1

   oplot,pcts[done],ys[done],psym=1,symsize=4.0

endelse


end

