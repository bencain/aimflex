FUNCTION set_union, a, b, count=count

; Note that the empty set is denoted by the value -1
;
; Code credit to:
; http://www.dfanning.com/tips/set_operations.html

  IF a[0] LT 0 THEN BEGIN
     count=n_elements(b)
     RETURN, b                                        ;A union NULL = a
  ENDIF
  IF b[0] LT 0 THEN BEGIN
     count=n_elements(a)
     RETURN, a                                        ;B union NULL = b
  ENDIF



  union=Where(Histogram([a,b],OMin = omin)) + omin
                                ; Return combined set

  if (size(a,/type) le 3) and (size(b,/type) le 3) then union=ceil(union)

  count=n_elements(union)
  return,union

END
