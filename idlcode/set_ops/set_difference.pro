FUNCTION set_difference, a, b  , count=count
; = a and (not b) = elements in A but not in B 
; Note that the empty set is denoted by the value -1
;
; Code credit to:
; http://www.dfanning.com/tips/set_operations.html

  if a[0] eq -1 then begin
     count=0
     return,-1
  endif

  if b[0] eq -1 then begin
     count=n_elements(a)
     return,a
  endif
  
  if size(a,/dimensions) eq 0 then a2=[a] else a2=a
  if size(b,/dimensions) eq 0 then b2=[b] else b2=b

  if ((n_elements(a2) eq 1) and (n_elements(b2) eq 1)) then begin
     if a2 eq b2 then begin
        count=0
        return,-1 
     endif else begin
        count=1
        return,a2
     endelse
  endif

  mina = Min(a2, Max=maxa)
  minb = Min(b2, Max=maxb)
  IF (minb GT maxa) OR (maxb LT mina) THEN BEGIN
     count=n_elements(a2)
     RETURN, a2                 ;No intersection...
  ENDIF

  r = Where((Histogram(a2, Min=mina, Max=maxa) NE 0) AND $
            (Histogram(b2, Min=mina, Max=maxa) EQ 0), count)
  IF count eq 0 THEN RETURN, -1 ELSE RETURN, r + mina
END
