pro disp_truefit, truecatfile, fitcatfile, skip=skip, start=start

truecat=read_data(truecatfile)
fitcat=read_data(fitcatfile)

ncat=min([n_elements(truecat[*,0]),n_elements(fitcat[*,0])])

if not keyword_set(start) then start=0L
if not keyword_set(skip) then skip=0L

dispok=1L
dispnum=start
psf=mk_gaussian_psf(3d)

while (dispok and (dispnum lt ncat)) do begin
   
   tp=truecat[dispnum,0:11]
   fp=fitcat[dispnum,0:11]

   parameter_status,tp-fp,dpnum,pnames=['dI0','dXc','dYc','dalpha','dEplus',$
                                        'dEcross','dg1','dg2','dG11','dG12',$
                                        'dG31','dG32'],$
                    label='Parameter Difference'
   parameter_status,fp,pnames=['I0','Xc','Yc','alpha','Eplus',$
                               'Ecross','g1','g2','G11','G12',$
                               'G31','G32'],label='Fit Parameters',fpnum

   chisq=fitcat[dispnum,12]
   dof=fitcat[dispnum,13]

   win=mk_window(6d*max([fp[3],tp[3]]))

   disp_scaled_image,mk_image(tp,85d,win,psf)-85d*win,tinum,title='True Image'
   disp_scaled_image,mk_image(fp,85d,win,psf)-85d*win,finum,title='Fit Image'
   xyouts,50d,50d,'chisq/dof = '+$
          strcompress(string(chisq)+'/'+string(dof),/remove_all),$
          /device,charsize=2d,charthick=2d
   xyouts,50d,25d,'Object = '+$
          strcompress(string(dispnum+1)+'/'+string(ncat),/remove_all),$
          /device,charsize=2d,charthick=2d

   clear
   dispok=get_num_resp('Continue? (1=yes, 0=no)', ok_resp=[1L,0L],default=1L)
   dispnum+=1+skip
endwhile

clear_windows

end
