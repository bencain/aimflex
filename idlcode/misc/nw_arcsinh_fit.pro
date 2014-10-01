;+
;NAME:
;  nw_arcsinh_fit
;PURPOSE:
;  scales the fits image by a degree of nonlinearity specified by user
;INPUTS:
;  colors      - (NXxNYx3) array that contains the R/G/B images
;OPTIONAL INPUTS:
;  nonlinearity- 'b'
;              - b=0 for linear fit, b=Inf for logarithmic
;              - default is 3
;KEYWORDS:
;
;OUTPUTS:
;  The image
;BUGS:
;  
;REVISION HISTORY:
;  11/07/03 written - wherry
;  11/12/03 changed radius - wherry
;-
FUNCTION nw_arcsinh_fit,colors,nonlinearity=nonlinearity
;set default nonlinearity
IF NOT n_elements(nonlinearity) THEN nonlinearity=3

R = colors[*,*,0]
G = colors[*,*,1]
B = colors[*,*,2]

dim = size(R,/dimensions)
NX = LONG(dim[0])
NY = LONG(dim[1])

radius = R+G+B
radius = radius+(radius eq 0)
IF (nonlinearity eq 0.) THEN BEGIN 
    val = radius
ENDIF ELSE BEGIN
    val = asinh(radius*nonlinearity)/nonlinearity
ENDELSE

fitted_colors = fltarr(NX,NY,3)
fitted_colors[*,*,0] = (R*val)/radius
fitted_colors[*,*,1] = (G*val)/radius
fitted_colors[*,*,2] = (B*val)/radius
RETURN,fitted_colors
END
