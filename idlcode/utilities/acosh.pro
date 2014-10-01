function acosh, x, neg_branch=neg_branch

; Inverse hyperbolic cosine

; cosh(y)^2-sinh(y)^2=1d
; sinh(y)=sqrt(cosh(y)^2-1)
;
; This fails on complex input.

on_error,2

y=asinh(sqrt(x^2-1d))

; acosh is double valued - allow for the negative branch to be
;                          returned also.
if keyword_set(neg_branch) then return,-y else return,y

end
