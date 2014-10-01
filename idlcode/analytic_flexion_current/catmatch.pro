FUNCTION CATMATCH, POINTS, CENTERS, RADII, COUNT=COUNT
  
; This function returns a list of elements of CENTERS/RADII which have
; an element of POINTS within RADII of CENTERS.  CENTERS is an Nx2
; element array, RADII is an N element array, and POINTS is Mx2.

npts=(size(points,/dimensions))[0]
nctr=min([(size(points,/dimensions))[0],n_elements(radii)])

flagged=-1
count=0

for i=0,npts-1 do begin
   ratio=sqrt((points[i,0]-centers[*,0])^2 + $
              (points[i,1]-centers[*,1])^2)/radii

   flagged=set_union(flagged,where(ratio le 1d),count=count)
endfor

return,flagged
end
