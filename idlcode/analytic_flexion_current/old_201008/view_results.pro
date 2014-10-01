pro view_results,tag

catalog=read_data(tag+'_fit.dat')

ncat=(size(catalog,/dimensions))[0]

minx=min(catalog[*,1])
miny=min(catalog[*,2])
maxx=max(catalog[*,1])
maxy=max(catalog[*,2])

midx=0.5d*(minx+maxx)
midy=0.5d*(miny+maxy)
widx=maxx-minx
widy=maxy-miny

xr=midx+widx*[-0.6d,0.6d]
yr=midy+widy*[-0.6d,0.6d]

for i=0,ncat-1 do begin

   objnum=long(catalog[i,0])

   if objnum+1 lt 10 then fill='000' else $
      if objnum+1 lt 100 then fill='00' else $
         if objnum+1 lt 1000 then fill='0' else fill=''

; Make the file names
   dataname=tag+'_data_'+fill+$
            strcompress(string(objnum),/remove_all)+'.fits'
   fitname=tag+'_fit_'+fill+$
           strcompress(string(objnum),/remove_all)+'.fits'
   residname=tag+'_resid_'+fill+$
             strcompress(string(objnum),/remove_all)+'.fits'


; Display the images, each in its own window
   dimg=mrdfits(dataname)
   fimg=mrdfits(fitname)
   rimg=mrdfits(residname)

   pos=where(dimg gt 0d)
   dimg[pos]-=catalog[i,19]
   fimg[pos]-=catalog[i,19]

   disp_scaled_image,dimg,dnum,title='Data Image',winsize=250
   disp_scaled_image,fimg,fnum,title='Fit Image',winsize=250
   disp_scaled_image,rimg,rnum,title='Residual Image',winsize=250
   
   loadct,13
   xyouts,0.1d,0.3d,'Mean = '+$
          strcompress(string(mean(rimg[pos])),/remove_all),/norm,charsize=2
   xyouts,0.1d,0.2d,'Median = '+$
          strcompress(string(median(rimg[pos])),/remove_all),/norm,charsize=2
   xyouts,0.1d,0.1d,'StdDev = '+$
          strcompress(string(stddev(rimg[pos])),/remove_all),/norm,charsize=2
   loadct,0


; Show the location of the galaxy
   if n_elements(p_num) eq 0 then begin
      window,xsize=250,ysize=250,/free,title='Galaxy Position'
      p_num=!d.window
   endif else wset,p_num
   plot,catalog[i,1],catalog[i,2],psym=2,xrange=xr,xstyle=1,$
        yrange=yr,ystyle=1,symsize=3
   loadct,13
   oplot,1979.281d,2122.3345d,psym=1,symsize=3
   loadct,0


; Show the lensing fields
   if n_elements(g_num) eq 0 then begin
      window,xsize=250,ysize=250,/free,title='g Field'
      g_num=!d.window
   endif else wset,g_num
   plot_field,catalog[i,*],g_num,xy_cols=[1,2],field_cols=[9,10],spin=2,$
              /standard,/std_label

   if n_elements(G1num) eq 0 then begin
      window,xsize=250,ysize=250,/free,title='G1 Field'
      G1num=!d.window
   endif else wset,G1num
   plot_field,catalog[i,*],G1num,xy_cols=[1,2],field_cols=[11,12],spin=1,$
              symscale=0.5d,/standard,/std_label
   
   if n_elements(G3num) eq 0 then begin
      window,xsize=250,ysize=250,/free,title='G3 Field'
      G3num=!d.window
   endif else wset,G3num
   plot_field,catalog[i,*],G3num,xy_cols=[1,2],field_cols=[13,14],spin=3,$
              symscale=0.5d,/standard,/std_label

   print,'*************************'
   print,'Object Number = '+strcompress(string(objnum),/remove_all)+' ('+$
         strcompress(string(i+1)+'/'+string(ncat),/remove_all)+')'
   print,'I0     = '+strcompress(string(catalog[i,3]),/remove_all)
   print,'Xc     = '+strcompress(string(catalog[i,4]),/remove_all)
   print,'Yc     = '+strcompress(string(catalog[i,5]),/remove_all)
   print,'alpha  = '+strcompress(string(catalog[i,6]),/remove_all)
   print,'eplus  = '+strcompress(string(catalog[i,7]),/remove_all)
   print,'ecross = '+strcompress(string(catalog[i,8]),/remove_all)
   print,'g1     = '+strcompress(string(catalog[i,9]),/remove_all)
   print,'g2     = '+strcompress(string(catalog[i,10]),/remove_all)
   print,'G11    = '+strcompress(string(catalog[i,11]),/remove_all)
   print,'G12    = '+strcompress(string(catalog[i,12]),/remove_all)
   print,'G31    = '+strcompress(string(catalog[i,13]),/remove_all)
   print,'G32    = '+strcompress(string(catalog[i,14]),/remove_all)
   print,'chisq  = '+strcompress(string(catalog[i,15]),/remove_all)
   print,'dof    = '+strcompress(string(catalog[i,16]),/remove_all)
   print,'n_iter = '+strcompress(string(catalog[i,18]),/remove_all)
   print,'*************************'

   ok=get_num_resp('Move to the next, or quit? (1=next, 0=quit)',default=1d)
   
   if ok eq 0d then break

endfor

clear_windows

end
