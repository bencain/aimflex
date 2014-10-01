pro fix_fits_nans,infile,outfile,value=value, nbr_ave=nbr_ave

  if n_elements(value) eq 0 then value=0d         ; Replacement value
  if n_elements(outfile) eq 0 then outfile=infile ; Default to overwrite if no outfile is given

  img=mrdfits(infile,0,hdr) ; Read in the image
  npix=n_elements(img)

  badpix=where(finite(img,/nan),nbad)
  
  if nbad gt 0 then begin
   
     img[badpix]=value          ; At least set them to zero
          
     if keyword_set(nbr_ave) then begin
        p1=badpix+1             ; Pick out the two neighbors
        p2=badpix-1
      
        kill1=where(p1 ge npix,nk1)
        kill2=where(p2 lt 0,nk2)
        if nk1 gt 0 then p1[kill1]=npix-1 ; Just use the last one
        if nk2 gt 0 then p2[kill2]=0      ; Just use the first one

        img[badpix]=(img[p1]+img[p2])*0.5d ; Average the neighbors
     endif

     spawn,'rm -f '+outfile
     mwrfits,img,outfile,hdr

  endif else print, 'No bad pixels, nothing changed'

end
   
