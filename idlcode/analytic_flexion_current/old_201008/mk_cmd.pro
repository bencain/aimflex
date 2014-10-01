pro mk_cmd, mag_arr, col_arr, mmin=mmin, mmax=mmax, cmin=cmin, cmax=cmax,$
            spread=spread, rs_frac=rs_frac, ngal=ngal

if not keyword_set(mmin) then mmin=0d
if not keyword_set(mmax) then mmax=6d
if not keyword_set(cmin) then cmin=-4d
if not keyword_set(cmax) then cmax=4d
if not keyword_set(spread) then spread=1d
if not keyword_set(rs_frac) then rs_frac=0.8d
if not keyword_set(ngal) then ngal=200L

nrs=long(ngal/3d)
nbulk=long(2*ngal/3d)

rs_mag=randomu(seed,nrs)*(mmax-mmin) + mmin
rs_col=((cmax-cmin)/(mmax-mmin))*(rs_mag-mmin)+spread*randomn(seed,nrs)+cmin

bulk_mag=mmax*(1-rs_frac)*randomu(seed,nbulk)+rs_frac*mmax
bulk_col=(cmax+cmin)*(0.5d)+(cmax-cmin)*(randomn(seed,nbulk)-0.5d)

mag_arr=[rs_mag,bulk_mag]
col_arr=[rs_col,bulk_col]

end
