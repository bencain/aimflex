;+
;NAME:
;  nw_scale_rgb
;PURPOSE:
;  mulitiplies the RGB image by their respective scales
;CALLING SEQUENCE:
;  nw_scale_rgb, colors, [scales=]
;INPUTS:
;  colors      - (NXxNYx3) array containing the R, G, and B
;OPTIONAL INPUTS:
;  scales      - (3x1) array to scale the R/G/B
;              - defaults are [4.9,5.7,7.8]
;KEYWORDS:
;  none
;OUTPUTS:
;  The RGB image 
;BUGS:
;  
;DEPENDENCIES:
;
;REVISION HISTORY:
;  11/07/03 written - wherry
;-
FUNCTION nw_scale_rgb,colors,scales=scales

;set default scales
IF NOT keyword_set(scales) THEN scales = [4.9,5.7,7.8]
;get dimensions
tmp= size(colors,/dimensions)
NX= LONG(tmp[0])
NY= LONG(tmp[1])

scaled_colors = fltarr(NX,NY,3)
FOR k=0,2 DO scaled_colors[*,*,k] = colors[*,*,k]*scales[k]

RETURN,scaled_colors
END
