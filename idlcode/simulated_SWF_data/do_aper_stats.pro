pro do_aper_stats,infile,outfile

  outexists=file_test(outfile)
  if outexists then openu,out,outfile,/get_lun,/append,width=200 else openw,out,outfile,/get_lun,width=200

; Read in the fluxes
  readcol,infile,ap,f,format='X,L,X,X,X,D'

  printf,out,ap[0],mean(f),stddev(f)

  close,out
  free_lun,out

end
