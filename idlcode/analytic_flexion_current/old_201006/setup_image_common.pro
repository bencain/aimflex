pro setup_image_common,npar,dims,sersic=sersic,$
                       mids=mids,scls=scls,pnames=pnames,$
                       clear=clear

common common_image,m,s,p

if keyword_set(clear) then begin
   delvarx,m
   delvarx,s
   delvarx,p
   
   delvarx,mids
   delvarx,scls
   delvarx,pnames
   
endif else begin

   maxs=dblarr(npar)
   mins=dblarr(npar)

   maxs[0:5]=[1d,$              ;I0
              0.25d*min(dims),$ ;Xc
              0.25d*min(dims),$ ;Yc
              0.5d*min(dims),$  ;alpha
              1d,$              ;E+
              1d]               ;Ex

   maxs[npar-6:npar-1]=[1d,$    ;g1
                        1d,$    ;g2
                        100d,$  ;G11
                        100d,$  ;G12
                        100d,$  ;G31
                        100d]   ;G32

   mins[0:5]=[1d-4,$             ;I0
              -0.25d*min(dims),$ ;Xc
              -0.25d*min(dims),$ ;Yc
              0.5d,$             ;alpha
              -1d,$              ;E+
              -1d]               ;Ex
   
   mins[npar-6:npar-1]=[-1d,$   ;g1
                        -1d,$   ;g2
                        -100d,$ ;G11
                        -100d,$ ;G12
                        -100d,$ ;G31
                        -100d]  ;G32

   pnames=strarr(npar)
   pnames[0:5]=['I0','X_c','Y_c','Alpha','E_+','E_x']
   pnames[npar-6:npar-1]=['g_1','g_2','G1_1','G1_2','G3_1','G3_2']
   
   if keyword_set(sersic) then begin
      maxs[6]=10d               ;index
      mins[6]=0.2d
      pnames[6]='n'
   endif

   mids=(maxs+mins)/2d
   scls=(maxs-mins)/2d
   
   m=mids
   s=scls
   p=pnames

endelse


end
