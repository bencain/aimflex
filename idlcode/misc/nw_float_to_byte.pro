;+
;NAME:
;  nw_float_to_byte
;PURPOSE:
;  Converts floats of an array to bytes
;INPUTS:
;  image       - image array
;OPTIONAL INPUTS:
;  none
;KEYWORDS:
;  none
;OUTPUTS:
;  The float-value image
;REVISION HISTORY:
;  10/03/03 written - wherry
;-
FUNCTION nw_float_to_byte,image
byte_image = byte((floor(image * 256.0) > 0) < 255)
RETURN, byte_image
END
