function cmd_selection, mag_arr, col_arr, nsigma=nsigma, redseq=redseq, $
                        removed=removed,pstag=pstag

; This function plots a color-magnitude diagram and takes user input
; to make cmd cuts on a cluster catalog.  This is intended to fit and
; remove the cluster red sequence.  We assume that the catalog has
; already been cleaned of stars with SIZE_SELECTION.PRO

if not keyword_set(nsigma) then nsigma=2d

; Select the galaxy objects which are ok to consider for RS fitting
objs=interactive_select(mag_arr,col_arr,aname='Magnitude',bname='Color',$
                        wtitle='Outlier Selection', not_sel=removed,$
                        prompt=$
                        'Select the usable objects, eliminating outliers.')

; Select out the galaxies for red sequence fitting
rsfit=interactive_select(mag_arr[objs],col_arr[objs],aname='Magnitude',$
                         bname='Color',wtitle='Select the Red Sequence',$
                         ainit_range=magrange, binit_range=colrange,$
                         prompt=$
                         'Select the objects for red sequence fitting.',$
                         psfile=pstag+'redseq.ps')

rsfit=objs[rsfit]

; Now fit the red sequence CMR.

cmr=poly_fit(mag_arr[rsfit],col_arr[rsfit],1,/double)
cmr_sig=stddev(col_arr[rsfit]-(cmr[0]+cmr[1]*mag_arr[rsfit]))

redseq=objs[where(abs(col_arr[objs]-$
                      (cmr[0]+cmr[1]*mag_arr[objs])) le nsigma*cmr_sig)]
not_redseq=objs[where(abs(col_arr[objs]-$
                          (cmr[0]+cmr[1]*mag_arr[objs])) gt nsigma*cmr_sig)]

; Plot the results
window,xsize=500,ysize=500,/free,title='Final selections'
plot,mag_arr[redseq],col_arr[redseq],psym=2,xrange=magrange,yrange=colrange,$
     xtitle='Magnitude',ytitle='Color',title='Magnitude vs Color'
oplot,mag_arr[not_redseq],col_arr[not_redseq],psym=1

; Plot the CMR with selection range
oplot,magrange,cmr[0]+magrange*cmr[1]
oplot,magrange,cmr[0]+magrange*cmr[1]+nsigma*cmr_sig,linestyle=2
oplot,magrange,cmr[0]+magrange*cmr[1]-nsigma*cmr_sig,linestyle=2

if keyword_set(pstag) then begin

   setup_ps,pstag+'_cmd_selection.ps'

   plot,mag_arr[redseq],col_arr[redseq],psym=2,xrange=magrange,$
        yrange=colrange,xtitle='Magnitude',ytitle='Color',$
        title='Magnitude vs Color'
   oplot,mag_arr[not_redseq],col_arr[not_redseq],psym=1
   
; Plot the CMR with selection range
   oplot,magrange,cmr[0]+magrange*cmr[1]
   oplot,magrange,cmr[0]+magrange*cmr[1]+nsigma*cmr_sig,linestyle=2
   oplot,magrange,cmr[0]+magrange*cmr[1]-nsigma*cmr_sig,linestyle=2
   
   setup_x

endif

return,not_redseq

end




