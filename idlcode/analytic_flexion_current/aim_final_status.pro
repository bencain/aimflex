pro aim_final_status, foms, iter_ind, peg_ind, fail_ind, maxiter=maxiter
if not keyword_set(maxiter) then maxiter=500

ftemp=where(foms[*,fail_ind] gt 0,nf)
ptemp=where(foms[*,peg_ind] gt 0,np)
itemp=where(foms[*,iter_ind] gt maxiter-1,ni)

fptemp=set_intersection(ftemp,ptemp,count=nfp)
fitemp=set_intersection(ftemp,itemp,count=nfi)
pitemp=set_intersection(ptemp,itemp,count=npi)
fpitemp=set_intersection(fptemp,itemp,count=nfpi)
anytemp=where((foms[*,fail_ind] gt 0) or $
              (foms[*,peg_ind] gt 0) or $
              (foms[*,iter_ind] gt maxiter-1),nany)

print,'Done!'
print,'Failed fits:         '+$
      strcompress(string(nf),/remove_all)
print,'Pegged fits:         '+$
      strcompress(string(np),/remove_all)
print,'Iterated out fits:   '+$
      strcompress(string(ni),/remove_all)
print,'Failed and pegged:   '+$
      strcompress(string(nfp),/remove_all)
print,'Pegged and iterated: '+$
      strcompress(string(npi),/remove_all)
print,'Failed and iterated: '+$
      strcompress(string(nfi),/remove_all)
print,'All three:           '+$
      strcompress(string(nfpi),/remove_all)
print,'Any of the three:    '+$
      strcompress(string(nany),/remove_all)
print,'Total # of fits:     '+$
      strcompress(string(n_elements(foms[*,0])),/remove_all)

end
