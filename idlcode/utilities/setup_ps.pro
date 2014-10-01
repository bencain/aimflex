pro setup_ps, psfile

; This is a shorthand to set the plotting to postscript.

  set_plot,'PS'
  device,filename=psfile,/color

end
