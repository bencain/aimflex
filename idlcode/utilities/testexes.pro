pro testexes


readcol,'~/himed.txt',lambda,orders,angles,coverage,slitlen,rnum

readcol,'~/hilo.txt',l,o,a,c,s,r
lambda = [lambda,l]
orders = [orders,o]
angles = [angles,a]
coverage = [coverage,c]
slitlen = [slitlen,s]
rnum = [rnum,r]
delvarx,l,o,a,c,s,r

readcol,'~/med.txt',l,o,a,c,s,r
lambda = [lambda,l]
orders = [orders,o]
angles = [angles,a]
coverage = [coverage,c]
slitlen = [slitlen,s]
rnum = [rnum,r]
delvarx,l,o,a,c,s,r

readcol,'~/lo.txt',l,o,a,c,s,r
lambda = [lambda,l]
orders = [orders,o]
angles = [angles,a]
coverage = [coverage,c]
slitlen = [slitlen,s]
rnum = [rnum,r]
delvarx,l,o,a,c,s,r


setup_ps,'~/orders_vs_lambda.ps'
plot,lambda,orders,/psym
setup_x

setup_ps,'~/angles_vs_lambda.ps'
plot,lambda,angles,/psym
setup_x

setup_ps,'~/speccov_vs_lambda.ps'
plot,lambda,coverage,/psym
setup_x

setup_ps,'~/slitlen_vs_lambda.ps'
plot,lambda,slitlen,/psym
setup_x

setup_ps,'~/rnum_vs_lambda.ps'
plot,lambda,rnum,/psym,/ylog
setup_x



end

