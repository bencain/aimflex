pro subcats, name, nobj_sub, ncats, seed=seed, tag=tag
  if not keyword_set(tag) then tag='sub'

  nobj_sub=long64(nobj_sub)
  ncats=long64(ncats)

  cat=read_data(name+'.dat',comment=com)
  nobj=n_elements(cat[*,0])

  if nobj le nobj_sub then stop,'Too many objects in the subcats - quitting.'
  if ncats lt 1 then stop,'Too few subcats requested'

 
; We'll subsample by source family.  The exact numbers won't be the
; same as nobj_sub but that's ok for now.
  
; find the total number of families
  nfam=long64(max(cat[*,21]))+1L

  i=0
  while i lt ncats do begin
     ntag=strcompress(string(i),/remove_all)
     subname=name+'_'+tag+'_'+ntag

     fam_use=long64(randomu(seed,nobj_sub,/double)*nfam)
     
; Kill any doubles and produce a sorted list of unique values.
     fam_use=fam_use[uniq(fam_use,sort(fam_use))]
     nfam_use=n_elements(fam_use)
     
; Where is the first image of each family
     starts=lon64arr(nfam_use)
     for j=0,nfam_use-1 do starts[j]=(where(cat[*,21] eq fam_use[j]))[0]
     
; How many images in each family?
     nimgs=cat[starts,22]
     
     if max(nimgs) gt 1 then begin

        subcat=dblarr(total(nimgs),24)
        substarts=total(nimgs,/cumul) - nimgs
     
        for j=0,nfam_use-1 do subcat[substarts[j]:substarts[j]+nimgs[j]-1,*]=cat[starts[j]:starts[j]+nimgs[j]-1,*]
        
        save_data,subcat,subname+'.dat',comment=com
        split_cats,subname,seed=seed

        i++

     endif
     
  endwhile

  print,'Done!'
end  
