function smc_selection, size_arr, mag_arr, col_arr, stars=stars, $
                        smc_outliers=smc_outliers, pstag=pstag

; This function plots color-size and magnitude-size diagrams and takes
; user input to make color/magnitude/size cuts on a galaxy catalog.

; Cut out the objects with excessively large sizes/bad colors.
ok=interactive_select(size_arr,col_arr,aname='Semimajor axis',$
                      bname='Color',not_sel=smc_outliers, prompt=$
                      'Select the objects with good size/color limits.',$
                      wtitle='Select Non-outliers')

; Get the star selection from size/color
cstars=interactive_select(size_arr[ok],col_arr[ok],aname='Semimajor axis',$
                          bname='Color',not_sel=cother,prompt=$
                          'Select the stars in size/color space.',$
                          wtitle='Select Stars',psfile=pstag+$
                          'sc_starsel.ps')

; Get the star selection from size/magnitude
mstars=interactive_select(size_arr[ok],mag_arr[ok],aname='Semimajor axis',$
                          bname='Magnitude',not_sel=mother,prompt=$
                          'Select the stars in size/magnitude space.',$
                          wtitle='Select Stars',psfile=pstag+$
                          'sm_starsel.ps')

; The STARS keyword is the intersection of the stars selected by color
; and by magnitude.
cstars=ok[cstars]
mstars=ok[mstars]
stars=long(set_union(mstars,cstars))

; The return value is anything that is not a star
if cother[0] ge 0 then cother=ok[cother]
if mother[0] ge 0 then mother=ok[mother]

other=long(set_intersection(cother,mother))


return,other

end


