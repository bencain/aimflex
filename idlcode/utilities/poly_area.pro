; $Id: poly_area.pro,v 1.6.4.1 1999/01/16 16:43:33 scottm Exp $
;
; Copyright (c) 1984-1999, Research Systems, Inc.  All rights reserved.
;       Unauthorized reproduction prohibited.

Function Poly_area,x,y, SIGNED=signed
;+
; NAME:
;	POLY_AREA
;
; PURPOSE:
;	Return the area of a polygon given the coordinates
;	of its vertices.
;
; CATEGORY:
;	Analytical Geometry
;
; CALLING SEQUENCE:
;	Result = POLY_AREA(X, Y)
;
; INPUTS:
;	It is assumed that the polygon has N vertices with N sides
;	and the edges connect the vertices in the order:
;
;	[(x1,y1), (x2,y2), ..., (xn,yn), (x1,y1)].
;
;	i.e. the last vertex is	connected to the first vertex.
;
;	X:	An N-element vector of X coordinate locations for the vertices.
;
;	Y:	An N-element vector of Y coordinate locations for the vertices.
;
; Keyword Inputs:
;	SIGNED = If set, returned a signed area. Polygons with edges
;	listed in counterclockwise order have a positive area, while those
;	traversed in the clockwise direction have a negative area.
; OUTPUTS:
;	POLY_AREA returns the area of the polygon.  This value is
;	positive, unless the SIGNED keyword is set and the polygon is
;	in clockwise order.
;
; COMMON BLOCKS:
;	None.
;
; SIDE EFFECTS:
;	None.
;
; RESTRICTIONS:
;	None.
;
; PROCEDURE:
;	The area is computed as:
;		Area = 	1/2 * [ x1y2 + x2y3 + x3y4 +...+x(n-1)yn + xny1 
;			- y1x2 - y2x3 -...-y(n-1)xn - ynx1)
;
; MODIFICATION HISTORY:
;	DMS, July, 1984.
;	DMS, Aug, 1996, Added SIGNED keyword.
;-
on_error,2                      ;Return to caller if an error occurs
n = n_elements(x)
if (n le 2) then message, 'Not enough vertices'
if n ne n_elements(y) then message,'X and Y arrays must have same size'

; Check type of arithmetic result, be sure we do things in doubles
if (size(x[0] + y[0]))[1] lt 4 then begin
    xx=double(x)		;Be sure its double
    a = total(xx*shift(y,-1) - y*shift(xx,-1))/2. ;This is it.
endif else a = total(x*shift(y,-1) - y*shift(x,-1))/2.  ;Already double

if keyword_set(signed) then return, a else return, abs(a)
end

