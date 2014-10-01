PRO CLEAN_IMG, IMAGEFILE, ERRIMGFILE, OBJIMGFILE, CATFILE1, CATFILE2, $
               OUTTAG=OUTTAG, SEED=SEED, $
               PRESELECTION=PRESELECTION,$
               FIT_OBJ=FIT_OBJ, RS_OBJ=RS_OBJ, $
               BAD_OBJ=BAD_OBJ, PT_OBJ=PT_OBJ,$
               N_FIT_OBJ=N_FIT_OBJ, N_RS_OBJ=N_RS_OBJ, $
               N_BAD_OBJ=N_BAD_OBJ, N_PT_OBJ=N_PT_OBJ,$
               NSIGMA=NSIGMA
if not keyword_set(outtag) then outtag='cleaned_' else outtag+='_'
if not keyword_set(nsigma) then nsigma=2d
; Read in the input
offset=1d

ext=0
img=mrdfits(imagefile,ext)
if n_elements(img) le 1 then begin
   ext++
   img=mrdfits(imagefile,ext)
endif
ext=0
err=mrdfits(errimgfile)
if n_elements(err) le 1 then begin
   ext++
   err=mrdfits(errimgfile,ext)
endif
ext=0
obj=mrdfits(objimgfile)
if n_elements(obj) le 1 then begin
   ext++
   obj=mrdfits(objimgfile,ext)
endif

dirtypix=dblarr(size(img,/dimensions))

outimg=img
outerr=err

if not keyword_set(preselection) then preselection=-1
toclean=-1
if get_num_resp('Do a CMS selection? (1=yes, 0=no)  ',default=1) then begin
   cat1=read_data(catfile1,size=catsize,/quiet)
   cat2=read_data(catfile2,/quiet)

   mag=cat2[*,8]
   col=cat1[*,8]-mag
   scl=cat2[*,4]

; Select which objects to clean
   fit_obj=aim_cms_select(col,mag,scl,$
                          n_fit_obj=n_fit_obj,$
                          rs_obj=rs_obj, n_rs_obj=n_rs_obj,$
                          pt_obj=pt_obj, n_pt_obj=n_pt_obj,$
                          badcm=bad_obj, n_badcm=n_bad_obj,$
                          nsigma=nsigma)

   toclean=set_union(rs_obj,pt_obj)
endif else begin
   cat1=read_data(catfile1,size=catsize,/quiet)
   cat2=read_data(catfile1,/quiet)
endelse

toclean=set_union(preselection,toclean,count=ntoclean)

; Clean the image

for i=0,ntoclean-1 do begin
   forloop_status,i,ntoclean,wnum,label='Cleaning objects...'

; Find the blob
   pos=floor(cat2[toclean[i],1:2])
   blob=blobmaker(obj,pos,min=0d,nblob=nblob)

; Mark all of the pixels to be cleaned. 
   if nblob gt 0 then dirtypix[blob]=1

endfor
forloop_status,0,0,wnum,/delete

; Clean the images
blob=where(dirtypix gt 0,nblob)
if nblob gt 0 then begin
   outimg[blob]-=obj[blob]
   outimg[blob]+=randomn(seed,nblob)*err[blob]
   outerr[blob]=sqrt(outerr[blob]^2+err[blob]^2)
endif

mwrfits,outimg,outtag+'img.fits'
mwrfits,outerr,outtag+'rms.fits'

END
