PRO PLOT_HIST, X, BINSIZE=BINSIZE, _Ref_Extra=extra, xrange=xrange,$
               yrange=yrange, norm=norm, overplot=overplot,xhatch=xhatch,$
               peak=peak,hatch_size=hatch_size,xlog=xlog

valuerange=max(x)-min(x)
if valuerange eq 0d then valuerange=0.1d

if not keyword_set(binsize) then binsize=valuerange/100d

ys=[0,0,histogram(x,binsize=binsize),0,0]

nbin=n_elements(ys)
xs=(dindgen(nbin)-2d)*binsize + min(x)

if keyword_set(xlog) then begin
   ok=where((xs gt 0d) and (ys gt 0d),nok)
   if nok gt 0 then begin
      xs=xs[ok]
      ys=ys[ok]
      nbin=nok
   endif
endif

if not keyword_set(xrange) then xrange=xs[[0,nbin-1]]

; Round the xrange to the nearest powers of ten that enclose the data.
if keyword_set(xlog) then begin
   xrange=(10d)^[floor(alog10(xrange[0])),ceil(alog10(xrange[1]))]

   xs=[max([xs[0]-binsize/2d,xrange[0]]),xs,xs[nbin-1]+binsize/2d]
   ys=[0,ys,0]
   binsize+=2

endif


if keyword_set(norm) then begin
   ys/=double(n_elements(x))
   if n_elements(yrange) ne 2 then yrange=[0d,min([1d,1.5d*max(ys)])]
endif else if n_elements(yrange) ne 2 then yrange=[min(ys),max(ys)*1.1d]

if keyword_set(overplot) then begin
   oplot,xs,ys,$
        psym=10,_Strict_Extra=extra
endif else begin
   plot,xs,ys,$
        psym=10,_Strict_Extra=extra,xrange=xrange,yrange=yrange,$
        xlog=keyword_set(xlog),xstyle=1,ystyle=1
endelse

if keyword_set(xhatch) then begin

   if not keyword_set(hatch_size) then begin
      hatch_size=max(ys)/double(max([nbin/2d,20]))
      
   endif else hatch_size=double(hatch_size)

   do_hor_hatch=where(ys gt hatch_size,n_do_hor)
   do_ver_hatch=where(ys gt 0d,n_do_ver)   
;   n_do_ver=0

   for i=0,n_do_hor-1 do begin
      n_hor_hatch=ceil(ys[do_hor_hatch[i]]/hatch_size)-1

      xlo=(xs[do_hor_hatch[i]-1]+xs[do_hor_hatch[i]])/2d
      xhi=(xs[do_hor_hatch[i]+1]+xs[do_hor_hatch[i]])/2d

      for j=0,n_hor_hatch-1 do begin

         oplot,[xlo,xhi],double([j,j]+1)*hatch_size

      endfor
   endfor

   for i=0,n_do_ver-1 do begin
      
      xlo=(xs[do_ver_hatch[i]-1]+xs[do_ver_hatch[i]])/2d
      xhi=(xs[do_ver_hatch[i]+1]+xs[do_ver_hatch[i]])/2d
      x0=(xs[do_ver_hatch[0]]+xs[do_ver_hatch[0]-1])/2d
      
      offset=(x0-xlo) mod hatch_size

      n_ver_hatch=floor((xhi-xlo-offset)/hatch_size)

      for j=0,n_ver_hatch-1 do begin

         oplot,(xlo+offset+double(j+1)*hatch_size)*[1d,1d],$
               [ys[do_ver_hatch[i]],0d]
      endfor
   endfor


endif


peak=max(ys)

END
