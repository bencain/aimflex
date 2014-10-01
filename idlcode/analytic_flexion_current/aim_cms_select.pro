function AIM_CMS_SELECT, COLOR, MAGNITUDE, SCALE, $
                         n_fit_obj=n_fit_obj, $
                         pt_obj=pt_obj, n_pt_obj=n_pt_obj, $
                         rs_obj=rs_obj, n_rs_obj=n_rs_obj, $
                         badcm=badcm, n_badcm=n_badcm,$
                         nsigma=nsigma
; This function uses the color, magnitude, and size arrays from the
; SExtractor catalogs (or elsewhere) and splits the objects into
; point sources, RS galaxies and fit objects.

if not keyword_set(nsigma) then nsigma=2d

; Make sure that the input is sane
nc=n_elements(color)
nm=n_elements(magnitude)
ns=n_elements(scale)

n_obj=min([nc,nm,ns])

if n_obj le 1 then begin
   fit_obj= -1
   pt_obj = -1
   rs_obj = -1
   badcm=-1
   n_fit_obj= 0
   n_pt_obj = 0
   n_rs_obj = 0
   n_badcm=0
   return,fit_obj
endif

c_arr=color[0:n_obj-1]
m_arr=magnitude[0:n_obj-1]
s_arr=scale[0:n_obj-1]

; An overall fitting loop
allok=(1 eq 0)
while not allok do begin

   clear_windows

; Check for bad colors/magnitudes
   cok=interactive_select(m_arr,c_arr,n_sel=n_cok,n_not_sel=n_badcm,not_sel=badcm,$
                          prompt='Select the region with well-defined color/magnitude.',$
                          aname='Magnitude',bname='Color',$
                          a_range=good_mrange,b_range=good_crange)
      
   clear_windows

; Now select on size/magnitude, if desired
   cut_pt_obj=(1 eq get_num_resp('Do you want to select out point sources? (1=yes 0=no):'))   
   if cut_pt_obj then begin
      ok=(1 eq 0)
      while (not ok) do begin
         pt_obj=interactive_select(s_arr[cok],m_arr[cok],$
                                   n_sel=n_pt_obj,not_sel=ext_obj,n_not_sel=n_ext_obj,$
                                   aname='Size',bname='Magnitude',$
                                   prompt='Select the point sources.',$
                                   a_range=pt_sizes,b_range=pt_mags)
         if (n_cok gt 0) and (n_pt_obj gt 0) then pt_obj=cok[pt_obj] else pt_obj=-1
         if (n_cok gt 0) and (n_ext_obj gt 0) then ext_obj=cok[ext_obj] else ext_obj=-1
         
         print,'By this selection: '
         print,strcompress(string(n_pt_obj)+'/'+string(n_cok),/remove_all)+$
               ' objects are point sources with ok colors/magnitudes.'
         print,strcompress(string(n_ext_obj)+'/'+string(n_cok),/remove_all)+$
               ' objects are extended sources with ok colors/magnitudes.'
         ok=get_num_resp('Is this okay? (1=yes, 0=no)')

         clear_windows
      endwhile
      
   endif else begin
      pt_obj=-1
      n_pt_obj=0
      ext_obj=cok
      n_ext_obj=n_cok
   endelse



; Check if we want to fit a red sequence

   plot,m_arr[ext_obj],c_arr[ext_obj],/psym,xtitle='Magnitude',ytitle='Color'

   fit_for_rs=(1 eq get_num_resp('Do you want to fit a red sequence? (1=yes 0=no):'))   
   if fit_for_rs then begin
      
; Select two points on the red sequence
   
      ok=(1 eq 0)
      while (not ok) do begin
   
         plot,m_arr[ext_obj],c_arr[ext_obj],/psym,xtitle='Magnitude',ytitle='Color'
         
         mmax=max(m_arr[ext_obj])
         mmin=min(m_arr[ext_obj])
         cmax=max(c_arr[ext_obj])
         cmin=min(c_arr[ext_obj])
         m1=mmin
         m2=mmax
         c1=cmin
         c2=cmax
         w=abs(c2-c1)
         print,'*******'
         print,'Enter two points, (m1,c1) and (m2,c2)'+$
               ' on the red sequence to start the RS selection:'
         read,m1,prompt='m1 = '
         read,c1,prompt='c1 = '
         oplot,[m1],[c1],color=fsc_color("Red"),psym=4,symsize=2,thick=3

         read,m2,prompt='m2 = '
         read,c2,prompt='c2 = '
         oplot,[m2],[c2],color=fsc_color("Red"),psym=4,symsize=2,thick=3

         read,w,prompt='Enter the initial guess for the width in color:  '
         oplot,[m1,m2,m2,m1,m1],[c1+w,c2+w,c2-w,c1-w,c1+w],color=1000,thick=2
         winit=w
         ok=get_num_resp('Is this okay? (1=yes, 0=no)')

         clear_windows
      endwhile



; Pick out the objects within w of the initial guess to do an
; iterative clipping fit of the RS.
      linear_rs=poly_fit([m1,m2],[c1,c2],1,/double)

; Do a sigma clipping loop
      ntest_prev=n_obj+10
      ok=(1 eq 0)
      while (not ok) do begin
         test=ext_obj[where(abs(c_arr[ext_obj] - $
                                (linear_rs[0]+linear_rs[1]*m_arr[ext_obj])) lt nsigma*w,ntest)]
         
         linear_rs=poly_fit(m_arr[test],c_arr[test],1,/double)
         w=stddev(c_arr[test] - (linear_rs[0]+linear_rs[1]*m_arr[test]))
         if ((ntest eq 0) or (ntest eq ntest_prev)) then break
         ntest_prev=ntest
      endwhile
      rs_obj=test
      n_rs_obj=ntest
   endif else begin
      rs_obj=-1
      n_rs_obj=0
   endelse

; Split out the fit objects.
   fit_obj=set_difference(lindgen(n_obj),$
                          set_union(rs_obj,set_union(pt_obj,badcm)),count=n_fit_obj)

   plot,m_arr[cok],c_arr[cok],/psym,xtitle='Magnitude',ytitle='Color'
   if n_rs_obj gt 0 then oplot,m_arr[rs_obj],c_arr[rs_obj],/psym,color=fsc_color("Red")
   if n_pt_obj gt 0 then oplot,m_arr[pt_obj],c_arr[pt_obj],/psym,color=fsc_color("Blue")
   print,'Fit',n_fit_obj
   print,'Pt ',n_pt_obj
   print,'RS ',n_rs_obj
   print,'bad',n_badcm

   allok=get_num_resp('Is this CMR selection okay? (1=yes, 0=no)')
endwhile

clear_windows
; Output to the screen a summary of the selection input
print,'**************************'
print,'**************************'
; Bad objects
print,'BADCM, N_BADCM'
print,'Objects with good color/magnitude were defined as:'
print,strcompress(string(good_crange[0]),/remove_all)+' < COLOR < '+$
      strcompress(string(good_crange[1]),/remove_all)
print,strcompress(string(good_mrange[0]),/remove_all)+' < MAGNITUDE < '+$
      strcompress(string(good_mrange[1]),/remove_all)

if cut_pt_obj then begin
   print,'**************************'
   print,'PT_OBJ, N_PT_OBJ'
   print,'Point sources were defined as:'
   print,strcompress(string(pt_sizes[0]),/remove_all)+' < ALPHA < '+$
         strcompress(string(pt_sizes[1]),/remove_all)
   print,strcompress(string(pt_mags[0]),/remove_all)+' < MAGNITUDE < '+$
         strcompress(string(pt_mags[1]),/remove_all)
endif

if fit_for_rs then begin
   print,'**************************'
   print,'RS_OBJ, N_RS_OBJ'
   print,'From the points ('+$
         strcompress(string(m1)+','+string(c1),/remove_all)+') and ('+$
         strcompress(string(m2)+','+string(c2),/remove_all)+'),' 
   print,'and an initial width of '+$
         strcompress(string(winit),/remove_all)+' the red sequence was found:'
   print,'SLOPE     = '+strcompress(string(linear_rs[1]),/remove_all)
   print,'INTERCEPT = '+strcompress(string(linear_rs[0]),/remove_all)
   print,'STD DEV   = '+strcompress(string(w),/remove_all)
   print,'Objects within 2*STDDEV of this RS were removed.'
endif
print,'**************************'
print,'FIT_OBJ, N_FIT_OBJ'
print,'All remaining objects are the fit objects.'
print,'**************************'
print,'N_FIT_OBJ  = '+strcompress(string(n_fit_obj),/remove_all)
print,'N_RS_OBJ   = '+strcompress(string(n_rs_obj),/remove_all)
print,'N_PT_OBJ   = '+strcompress(string(n_pt_obj),/remove_all)
print,'N_BADCM = '+strcompress(string(n_badcm),/remove_all)
print,'**************************'
print,'**************************'

ok=get_num_resp('Enter 1 when ready to continue.',ok_resp=1,default=1)

return,fit_obj

end


