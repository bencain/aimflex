function catalog_match, master, to_match, unmatched=unmatched,$
                        index_col=index_col, xy_cols=xy_cols, tol=tol

; This function matches the catalog index of objects in the master
; catalog to an unmatched catalog.  MASTER and TO_MATCH are assumed
; to be of the same column format, though they could have different
; numbers of objects.  For each object in TO_MATCH we find the
; closest object in MASTER, and if that object is within the given
; distance tolerance, TOL, the index of the master object is assigned
; to the unmatched object. If there is no master object within TOL of
; the unmatched object, then the index for that unmatched object is
; set to -1.  Those objects are returned in the UNMATCHED named
; variable. 

if not keyword_set(index_col) then index_col=0
if not keyword_set(xy_cols) then xy_cols=[1,2]
if not keyword_set(tol) then tol=2d

master_ind=master[*,index_col]
master_x=master[*,xy_cols[0]]
master_y=master[*,xy_cols[1]]

matched=to_match

n_to_match=(size(to_match,/dimensions))[0]

for i=0,n_to_match-1 do begin
   
   dists=sqrt((to_match[i,xy_cols[0]]-master_x)^2+$
              (to_match[i,xy_cols[1]]-master_y)^2)

   best=where(dists eq min(dists))

   if dists[best[0]] lt tol then begin
      matched[i,index_col]=master[best[0],index_col] 
   endif else begin
      matched[i,index_col]=-1
   endelse

endfor

neg=where(matched[*,index_col] lt 0,nneg,complement=pos,ncomplement=npos)

if nneg gt 0 then unmatched=matched[neg,*] else unmatched=-1
matched=matched[pos,*]

order=sort(matched[*,index_col])

return,matched[order,*]

end
