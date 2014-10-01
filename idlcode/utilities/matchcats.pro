pro matchcats

path='/Users/bcain/analysis_code/dropouts/F850LP_dropouts_MACS0744/'

files=path+['F435W','F606W','F775W','F850LP','F105W','F110W','F125W','F140W','F160W']+'.cat'

ncat=n_elements(files)

readcol,map,mape,mau,maue,x,y,ra,dec

nobj=n_elements(map)


MAG_APER=dblarr(nobj,ncat)
MAGERR_APER=dblarr(nobj,ncat)
MAG_AUTO=dblarr(nobj,ncat)
MAGERR_AUTO=dblarr(nobj,ncat)
X_IMAGE=dblarr(nobj,ncat)
Y_IMAGE=dblarr(nobj,ncat)
ALPHA_J2000=dblarr(nobj,ncat)
DELTA_J2000=dblarr(nobj,ncat)

mag_aper[*,0]=map
magerr_aper[*,0]=mape
mag_auto[*,0]=mau
magerr_auto[*,0]=maue
x_image[*,0]=x
y_image[*,0]=y
alpha_j2000[*,0]=ra
delta_j2000[*,0]=dec

for i=1,ncat-1 do begin
   readcol,map,mape,mau,maue,x,y,ra,dec

   mag_aper[*,i]=map
   magerr_aper[*,i]=mape
   mag_auto[*,i]=mau
   magerr_auto[*,i]=maue
   x_image[*,i]=x
   y_image[*,i]=y
   alpha_j2000[*,i]=ra
   delta_j2000[*,i]=dec
endfor

used=dblarr(nobj,ncat)


      


tol=2d
