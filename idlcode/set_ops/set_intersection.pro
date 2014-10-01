FUNCTION set_intersection, a, b, count=count

; Note that the empty set is denoted by the value -1
;
; Code credit to:
; http://www.dfanning.com/tips/set_operations.html
  if n_elements(a) eq 0 then begin
     count=0
     return,-1
  endif
  if n_elements(b) eq 0 then begin
     count=0
     return,-1
  endif

  if (n_elements(a) eq 1) or (n_elements(b) eq 1) then $
     return,set_intersection([a,a],[b,b],count=count)
  minab = Min(a, Max=maxa) > Min(b, Max=maxb) 
                                ;Only need intersection of ranges 
  maxab = maxa < maxb
  
                                ; If either set is empty, or
                                ; their ranges don't intersect:
                                ; result = NULL 
  count=0
  IF maxab LT minab OR maxab LT 0 THEN RETURN, -1
  r = Where((Histogram(a, Min=minab, Max=maxab) NE 0) AND  $
            (Histogram(b, Min=minab, Max=maxab) NE 0), count)
  
  IF count EQ 0 THEN RETURN, -1 ELSE RETURN, r + minab
END
