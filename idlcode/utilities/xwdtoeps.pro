PRO XWDTOEPS, filename

array = READ_XWD(filename, r, g, b); Read the XWD file. Pixel intensity information is stored in the variable `array'. Values to reconstruct the color table are stored in `r', `g', and `b'.

TVLCT, r,g,b	; Reconstruct the color table.

TV, array	; Display the image in an IDL window.

s = SIZE(array)	; Find the size of the picture. The width of the picture (in pixels) is stored in s[1]. The height of the picture is stored in s[2]:

fl = STRLEN(filename)	; Take the `xwd' (for X Windows Dump) extension off of the old filename and replace it with `eps'.

filename = STRMID(filename, 0, fl-4)

filename = filename + '.eps'

PRINT, 'Making file: ', filename

PRINT, s

SET_PLOT, 'ps'	; Set the plotting device to PostScript.

DEVICE, /ENCAPSUL, BITS_PER_PIXEL=8, /COLOR, $

   FILENAME=filename, XSIZE=S[1]/1000., $

   YSIZE=S[2]/1000.	; Use the DEVICE procedure to make the output encapsulated, 8 bits, color, and only as wide and high as it needs to be to contain the XWD image:

TV, array	; Write the image to the file.

DEVICE, /CLOSE	; Close the file.

SET_PLOT, 'x'	; Return plotting to X Windows.

END

