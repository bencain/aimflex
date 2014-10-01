FUNCTION CATALOG_SELECT, MAG, COLOR, SIZE, PSTAG=PSTAG, TOEPS=TOEPS,$
                         SEL_UPPERONLY=SEL_UPPERONLY, NSIGMA=NSIGMA,$
                         OBJ=OBJ, N_OBJ=N_OBJ,$     ; Objs detected in both
                         GALS=GALS, N_GALS=N_GALS,$ ; Extended sources
                         PTS=PTS, N_PTS=N_PTS,$     ; Point sources
                         REDSEQ=REDSEQ, N_REDSEQ=N_REDSEQ,$  
                                                    ; Red-Seq galaxies
                         NOTSEL=NOTSEL, N_NOTSEL=N_NOTSEL,$
                                                    ; Non-RS gals not 
                                                    ; selected
                         NONDET=NONDET, N_NONDET=N_NONDET,$
                                                    ; Objs with bad colors
                         N_GOOD=N_GOOD,$
                         MAGTAG=MAGTAG, COLORTAG=COLORTAG,$
                         LIMITS=LIMITS
                                ; All objects in OBJ which are not in
                                ; REDSEQ or NOTSEL will be the return
                                ; value. The indices are referenced to
                                ; the original MAG, COLOR, and SIZE
                                ; arrays. LIMITS will contain the
                                ; selection limits. 

; This is a new version of the catalog selection function.  I plan for
; this to be a little more user-friendly and robust.

if not keyword_set(pstag) then pstag='catalog_select'
if not keyword_set(nsigma) then nsigma=2d
if not keyword_set(magtag) then magtag='Magnitude'
if not keyword_set(colortag) then colortag='Color'

; Make sure that there are the same number of elements across the board.
n_total=min([n_elements(mag),$
            n_elements(color),$
            n_elements(size)])

mag=mag[0:n_total-1]
color=color[0:n_total-1]
size=size[0:n_total-1]

; Cut out the objects with huge colors, which are non-detections in
; the blue frame.

nondet=where(color gt 90.0d,n_nondet,complement=obj,ncomplement=n_obj)

;Select out the point sources by size/magnitude

pts=interactive_select(mag[obj],size[obj],aname='Magnitude',$
                       bname='sqrt(AB)',wtitle='Point Source Selection',$
                       not_sel=gals,n_not_sel=n_gals,n_sel=n_pts,$
                       prompt='Select the objects which are point sources.',$
                       a_range=pt_magrng,b_range=pt_sizerng)

pts=obj[pts]
gals=obj[gals]


limits=[mean(pt_magrng),pt_magrng[1]-pt_magrng[0],$
        mean(pt_sizerng),pt_sizerng[1]-pt_sizerng[0]]

; Fit and cut out the red sequence.
clear
pick_rs=get_num_resp(['Select a red sequence?',$
                          '  0 = no',$
                          '  1 = yes'],default=0)
magrange=[min(mag),max(mag)]
colrange=[min(color),max(color)]

if pick_rs then begin
   rsfit=interactive_select(mag[gals],color[gals],aname='Magnitude',$
                            bname='Color',wtitle='Select the Red Sequence',$
                            ainit_range=magrange, binit_range=colrange,$
                            prompt=$
                            'Select the objects for red sequence fitting.',$
                            a_range=rs_magrng,b_range=rs_colrng)
   limits=[limits,$
           mean(rs_magrng),rs_magrng[1]-rs_magrng[0],$
           mean(rs_magrng),rs_colrng[1]-rs_colrng[0]]


   rsfit=gals[rsfit]

   cmr=poly_fit(mag[rsfit],color[rsfit],1,/double)
   cmr_sig=stddev(color[rsfit]-(cmr[0]+cmr[1]*mag[rsfit]))

   cmr_color=cmr[0]+cmr[1]*mag[gals]


   redseq=where(abs(color[gals]-cmr_color) lt nsigma*cmr_sig,n_redseq,$
                complement=notrs,ncomplement=n_notrs)

   redseq=gals[redseq]
   notrs=gals[notrs]
endif else begin
   redseq=-1
   n_redseq=0
   notrs=gals
   n_notrs=n_gals
endelse

; Pick a maximum magnitude (if desired).
clear
print,'Current galaxy magnitude limit: '+$
      strcompress(string(max(mag[gals])),/remove_all)
pick_maxmag=get_num_resp(['Select a new magnitude limit?',$
                          '  0 = no',$
                          '  1 = yes'],default=0)
if pick_maxmag eq 0 then maxmag=max(mag)+10d
if pick_maxmag eq 1 then begin
   read,maxmag,prompt='Enter magnitude limit:  '
   limits=[limits,maxmag]
endif

if keyword_set(sel_upperonly) then begin
   good=where((color[notrs] gt cmr_color[notrs]) and (mag[notrs] le maxmag),$
                                                      n_good,$
                                                      complement=notsel,$
                                                      ncomplement=n_notsel)
   good=notrs[good]
   notsel=notrs[notsel]
endif else begin
   good=where(mag[notrs] le maxmag,n_good,complement=notsel,$
              ncomplement=n_notsel)
   if n_good gt 0 then good=notrs[good]
   if n_notsel gt 0 then notsel=notrs[notsel]
endelse

; Plot the results

window,xsize=500,ysize=500,/free,title='Final selections'

; Plot the RS, the good galaxies, the unselected galaxies and the
; point sources.
plot,mag[good],color[good],psym=1,xrange=magrange,yrange=colrange,$
     xstyle=1,xtitle='Magnitude',ytitle='Color',title='Magnitude vs Color'
if n_redseq gt 0 then oplot,mag[redseq],color[redseq],psym=7
if n_notsel gt 0 then oplot,mag[notsel],color[notsel],psym=5
if n_pts gt 0 then oplot,mag[pts],color[pts],psym=2

; Plot the CMR with selection range and magnitude limit
if n_redseq gt 0 then begin
   oplot,magrange,cmr[0]+magrange*cmr[1]
   oplot,magrange,cmr[0]+magrange*cmr[1]+nsigma*cmr_sig,linestyle=2
   oplot,magrange,cmr[0]+magrange*cmr[1]-nsigma*cmr_sig,linestyle=2
   oplot,[maxmag,maxmag],colrange,linestyle=2
endif

if keyword_set(pstag) then begin
   
   setup_ps,pstag+'_cmd_sel.ps'

; Replot for the PS file
   plot,mag[good],color[good],psym=1,xrange=magrange,yrange=colrange,$
        xstyle=1,xtitle=magtag,ytitle=colortag,title='',$
        xmargin=[30,30], ymargin=[15,15],$
        charsize=3.0,charthick=3.0,thick=3.0,xthick=3.0,ythick=3.0
   if n_redseq gt 0 then oplot,mag[redseq],color[redseq],psym=7
   if n_notsel gt 0 then oplot,mag[notsel],color[notsel],psym=5
   if n_pts gt 0 then oplot,mag[pts],color[pts],psym=2

; Plot the CMR with selection range and magnitude limit
   if n_redseq gt 0 then begin
      oplot,magrange,cmr[0]+magrange*cmr[1]
      oplot,magrange,cmr[0]+magrange*cmr[1]+nsigma*cmr_sig,linestyle=2
      oplot,magrange,cmr[0]+magrange*cmr[1]-nsigma*cmr_sig,linestyle=2
      oplot,[maxmag,maxmag],colrange,linestyle=2
   endif

   setup_x

   if keyword_set(toeps) then ps2eps,pstag+'_cmd_sel.ps',/delete

endif


return,good
end
