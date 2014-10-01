FUNCTION SIGMA_CLIP,ARR, NSIGMA, NCLIP_MAX=NCLIP_MAX, NCLIP_USED=NCLIP_USED,$
                    MEDIAN=MEDIAN, NOK=NOK

; This function does a recursive sigma clipping of the array ARR.  The
; elements of the array further than NSIGMA standard deviations from
; the mean (or the median, if the MEDIAN keyword is set) will be
; removed and the mean/standard deviation will be recalculated until
; there are either no elements left to remove (all are within NSIGMA
; standard deviations) or until NCLIP clippings have occurred.

if not keyword_set(nclip_max) then nclip_max=5

nok=n_elements(arr)
ok=lindgen(nok)

for i=0, nclip_max-1 do begin

   nclip_used=i+1

   if keyword_set(median) then ctr=median(arr[ok]) else ctr=mean(arr[ok])
   sig=stddev(arr[ok])

   newok=where(abs(arr[ok]-ctr) lt nsigma*sig,nnewok)
   if nnewok gt 0 then newok=ok[newok]
   
   if (nnewok eq nok) or (nnewok eq 0) then begin
      ok=newok
      nok=nnewok
      return,ok
   endif else begin
      ok=newok
      nok=nnewok
   endelse
endfor

return,ok

end

   

   
