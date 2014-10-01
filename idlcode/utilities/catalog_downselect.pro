pro catalog_downselect, infiles, keepers, outfiles=outfiles, outtag=outtag

nin=n_elements(infiles)

if not keyword_set(outtag) then outtag='keep'

if n_elements(outfiles) eq 0 then begin
   outfiles=strarr(nin)
   for i=0,nin-1 do begin
      split=strsplit(infiles[i],'.',/extract,count=nsplit)
      outfiles[i]+=split[0]
      for j=1,nsplit-2 do outfiles[i]+='.'+split[j]
      outfiles[i]+='_'+outtag
      if nsplit gt 1 then outfiles[i]+='.'+split[nsplit-1]
   endfor
endif

nout=min([nin,n_elements(outfiles)])

for i=0,nout-1 do begin
   cat=read_data(infiles[i],comment=com,/quiet)
   save_data,cat[keepers,*],outfiles[i],comment=com
   print,infiles[i]+' ===> '+outfiles[i]
endfor

end
