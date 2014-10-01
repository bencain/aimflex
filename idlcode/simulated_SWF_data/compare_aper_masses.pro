pro compare_aper_masses


path='/Users/bcain/analysis_code/swunited_rewrite/simdata/simulation-2014-01-13/simulation01/'

readcol,path+'masses.txt',file,region,x,y,ra,dec,rpix,rwcs,m,format='A,L,D,D,D,D,D,D,D'

truth=where(strmatch(file,'*kappa.fits'))
recon=where(strmatch(file,'*kappa_*.fits'))


ft=  file[truth]
regt=region[truth]
xt=  x[truth]
yt=  y[truth]
rat= ra[truth]
dect=dec[truth]
rpt= rpix[truth]
rwt= rwcs[truth]
mt=  m[truth]

fr=  file[recon]
regr=region[recon]
xr=  x[recon]
yr=  y[recon]
rar= ra[recon]
decr=dec[recon]
rpr= rpix[recon]
rwr= rwcs[recon]
mr=  m[recon]

ok=where(mt gt 0d)


plot,mt[ok],mr[ok]/mt[ok],/psym

best=ok[where(abs(mt[ok]-mr[ok])/mt[ok] lt 1d)]

print,fr[best],regr[best]


zeroes=where(region eq 0)  
lasts=[zeroes[1:*]-1,n_elements(ra)-1]
n=region[lasts]+1


end
