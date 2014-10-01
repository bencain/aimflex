function sign, X, DOUBLE=DOUBLE
  
; This function returns an array which is the same size as X which
; takes on the value 1 where X[i] >= 0 and -1 where X[i] < 0. If the
; DOUBLE keyword is set, then that is the type returned.  Otherwise
; the return value is of long type.

sign=(2L)*(X ge 0)-1L

if keyword_set(double) then sign=double(sign)

return,sign

end
