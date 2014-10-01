function cat_overlap, xya1, xya2, factor=factor, noverlap=noverlap

; This procedure finds lists of which elements from array xya1
; overlap xya2 by factor times the radius in the (i,2) element of
; either.

if not keyword_set(factor) then factor=2d
single=0

; One array calling
if n_elements(xya2) eq 0 then begin
   xya2=xya1
   single=1
endif

s1=size(xya1,/dimensions)
s2=size(xya2,/dimensions)

if (n_elements(s1) ne 2) or (n_elements(s1) ne 2) then begin
   noverlap=0
   return,-1
endif

if (s1[1] ne 3) or (s2[1] ne 3) then begin
   noverlap=0
   return,-1
endif

bad=-1L

x1=xya1[*,0]
y1=xya1[*,1]
a1=xya1[*,2]

x2=xya2[*,0]
y2=xya2[*,1]
a2=xya2[*,2]

; See which objects overlap.  We only return things on array 1 which
; are bad. So we want to loop through each of the elements on array 2
; and see which elements of array 1 are overlapped.

for i=0L,s2[0]-1 do begin
   r=sqrt((x1-x2[i])^2+(y1-y2[i])^2)
   lim=a1

   switchlim=where(lim lt a2[i])
   if switchlim[0] ge 0 then lim[switchlim]=a2[i]
   lim*=factor

; If there is only one input array, cut out the identity match
   use=lindgen(s1[0])
   if single then use=set_difference(use,i)

   newbad=where(r[use] lt lim)
   if newbad[0] ge 0 then newbad=use[newbad]

 ;  print,newbad
   bad=set_union(bad,newbad)
 ;  print,bad
endfor


if bad[0] lt 0 then noverlap=0 else noverlap=n_elements(bad)
return,bad

end
