FUNCTION POLY_FIT, X, Y, NDEGREE, YFIT, YBAND, SIGMA, CORRM, DOUBLE=double
;+
; NAME:
;	POLY_FIT
;
; PURPOSE:
;	Perform a least-square polynomial fit with optional error estimates.
;
;	This routine uses matrix inversion.  A newer version of this routine,
;	SVDFIT, uses Singular Value Decomposition.  The SVD technique is more
;	flexible, but slower.
;
;	Another version of this routine, POLYFITW, performs a weighted
;	least square fit.
;
; CATEGORY:
;	Curve fitting.
;
; CALLING SEQUENCE:
;	Result = POLY_FIT(X, Y, NDegree [,Yfit, Yband, Sigma, CORRM] )
;
; INPUTS:
;	X:	The independent variable vector.
;
;	Y:	The dependent variable vector, should be same length as x.
;
;     NDegree:	The degree of the polynomial to fit.
;
; OUTPUTS:
;	POLY_FIT returns a vector of coefficients with a length of NDegree+1.
;
; OPTIONAL OUTPUT PARAMETERS:
;	Yfit:	The vector of calculated Y's.  These values have an error 
;		of + or - Yband.
;
;	Yband:	Error estimate for each point = 1 sigma
;
;	Sigma:	The standard deviations of the returned coefficients.
;
;	Corrm:	Correlation matrix of the coefficients.
;
; Keyword Parameters:
;	DOUBLE = if set, force computations to be in double precision.
; COMMON BLOCKS:
;	None.
;
; SIDE EFFECTS:
;	None.
;
; MODIFICATION HISTORY:
;	Written by: George Lawrence, LASP, University of Colorado,
;		December, 1981.
;
;	Adapted to VAX IDL by: David Stern, Jan, 1982.
;       Modified:    GGS, RSI, March 1996
;                    Corrected a condition which explicitly converted all
;                    internal variables to single-precision float.
;                    Added support for double-precision inputs.
;                    Added a check for singular array inversion.
;		     SVP, RSI, June 1996
;                     Changed A to Corrm to match IDL5.0 docs.
;                    S. Lett, RSI, December 1997
;                     Changed inversion status check to check only for
;                     numerically singular matrix.
;                    S. Lett, RSI, March 1998
;                     Initialize local copy of the independent variable
;                     to be of type DOUBLE when working in double precision.
;
;-
	ON_ERROR,2		;RETURN TO CALLER IF ERROR

	N = N_ELEMENTS(X) 	;SIZE
	IF N NE N_ELEMENTS(Y) THEN $
	  message,'X and Y must have same # of elements'
;
	M = NDEGREE + 1L ;# OF ELEMENTS IN COEFF VEC.
;

  sx = size(x)
  sy = size(y)
  if n_elements(double) eq 0 then $	;True if working in double prec
	double = (sx[sx[0]+1] eq 5) or (sy[sy[0]+1] eq 5)

  if double then begin
	CORRM = DBLARR(M,M)		;Operands are double, COEFF MATRIX
	B = DBLARR(M)		;WILL CONTAIN SUM Y * X^J
	Z = Replicate(1.0d0, N)
        xx = double(x)
  endif else begin
        CORRM = fltarr(M,M)         ;COEFF MATRIX
        B = fltarr(M)           ;WILL CONTAIN SUM Y * X^J
        Z = Replicate(1.0, N)
        xx = float(x)
  endelse
;
	B[0] = TOTAL(Y, DOUBLE = double)
	CORRM[0,0] = N
;
	FOR P = 1,2*NDEGREE DO BEGIN ;POWER LOOP.
		Z=Z*XX			;Z IS NOW X^P
                ;B IS SUM Y*X\X^J
		IF P LT M THEN B[P] = TOTAL(Y*Z)
		SUM = TOTAL(Z)
		FOR J= 0 > (P-NDEGREE), NDEGREE < P DO CORRM[J,P-J] = SUM
	  END			;END OF P LOOP.
;
	CORRM = INVERT(CORRM, status)	;INVERT MATRIX.
        if status eq 1 then message, "Singular matrix detected."
;
;			IF CORRM IS MULTIPLIED BY SIGMA SQUARED, IT IS THE
;			CORRELATION MATRIX.
;
	C = b # CORRM	;Get coefficients
;
	IF (N_PARAMS(0) LE 3) THEN RETURN,C	;EXIT IF NO ERROR ESTIMATES.
;
	YFIT = REPLICATE(C[0], N)  ;Initial Yfit
	FOR K = 1,NDEGREE DO YFIT = YFIT + C[K]*(XX^K) ;FORM YFIT.
;
	IF (N_PARAMS(0) LE 4) THEN RETURN,C	;EXIT IF NO ERROR ESTIMATES.
;
	IF N GT M THEN $ ;COMPUTE SIGMA
		SIGMA = TOTAL((YFIT-Y)^2) / (N-M) ELSE	SIGMA = 0.
;
	CORRM=CORRM* SIGMA		;GET CORRELATION MATRIX
;
	SIGMA = SQRT(SIGMA)
	YBAND = FLTARR(N)+ CORRM[0,0]	;SQUARED ERROR ESTIMATES
;
	FOR P = 1,2*NDEGREE DO BEGIN
	  Z = XX ^ P
	  SUM = 0.
	  FOR J=0 > (P - NDEGREE), NDEGREE < P DO SUM = SUM + CORRM[J,P-J]
	  YBAND = YBAND + SUM * Z ;ADD IN ERRORS.
	END		;END OF P LOOP
	YBAND = SQRT(ABS(YBAND))	;ERROR ESTIMATES
	RETURN,C
END
