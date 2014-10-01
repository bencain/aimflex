pro test_trunc

mu_one=mrdfits('/Users/bcain/analysis_code/swunited_rewrite/simdata/real3/real3_mu.fits')
k_one=mrdfits('/Users/bcain/analysis_code/swunited_rewrite/simdata/real3/real3_kappa.fits')
g1_one=mrdfits('/Users/bcain/analysis_code/swunited_rewrite/simdata/real3/real3_gamma1.fits')
g2_one=mrdfits('/Users/bcain/analysis_code/swunited_rewrite/simdata/real3/real3_gamma2.fits')

mu_sub=mrdfits('/Users/bcain/analysis_code/swunited_rewrite/simdata/real2/real2_mu.fits')
k_sub=mrdfits('/Users/bcain/analysis_code/swunited_rewrite/simdata/real2/real2_kappa.fits')
g1_sub=mrdfits('/Users/bcain/analysis_code/swunited_rewrite/simdata/real2/real2_gamma1.fits')
g2_sub=mrdfits('/Users/bcain/analysis_code/swunited_rewrite/simdata/real2/real2_gamma2.fits')

z=dindgen(50)/5d + 0.31d 
wt=redshift_to_weight(z,0.3d) 

mu_sub_ave=wt*0d
mu_one_ave=wt*0d
mu_sig=wt*0d
mu_ratio_ave=wt*0
mu_ratio_sig=wt*0

gsq_one=g1_one^2+g2_one^2
gsq_sub=g1_sub^2+g2_sub^2


for i=0,49 do begin
   mstmp=1d/( (1d - wt[i]*k_sub)^2 - wt[i]^2*gsq_sub )
   motmp=1d/( (1d - wt[i]*k_one)^2 - wt[i]^2*gsq_one )

   bad1=where(mstmp gt 100,nb1)
   if nb1 gt 0 then mstmp[bad1]=100d
   bad2=where(mstmp lt -100,nb2)
   if nb2 gt 0 then mstmp[bad2]=-100d

   bad1=where(motmp gt 100,nb1)
   if nb1 gt 0 then motmp[bad1]=100d
   bad2=where(motmp lt -100,nb2)
   if nb2 gt 0 then motmp[bad2]=-100d


   mu_sub_ave[i]=mean(abs(mstmp))
   mu_one_ave[i]=mean(abs(motmp))

   mu_sig[i]=stddev(abs(mstmp)-abs(motmp))
   mu_ratio_ave[i]=mean(abs(mstmp/motmp))
   mu_ratio_sig[i]=stddev(abs(mstmp/motmp))
endfor


!p.multi=[0,3,2]

plot,z,mu_sub_ave
oplot,z,mu_one_ave

plot,z,mu_ratio_ave
plot,z,mu_ratio_sig

plot_hist,abs(mstmp)
plot_hist,abs(motmp)

;!p.multi=0
pts=floor(1d6*randomu(seed,300))
plot,abs(mstmp[pts]),abs(mstmp[pts]/motmp[pts]),/psym

print,correlate(abs(mstmp),abs(mstmp/motmp))
end

x=2d*randomu(seed,500,2,/double) - 1d
x[*,0]=(dindgen(500)/500d + 0.0001d)*20d
x[*,1]=0d
c=[0d,0d]
b=0.1d
s=0.05d
t=0.3d
q=1d
a=0d

p=[c,b,s,t,q,a]

L=nie_trunc_lens(x,p)

pp=p[[0,1,2,3,5,6]]
pn=p[[0,1,2,4,5,6]]
print,pp
print,pn
Lp=nie_lens(x,pp)
Ln=nie_lens(x,pn)

r=sqrt(total(x^2,2,/double))

help
F=sqrt(total(L[*,6:7]^2,2,/double))
Fp=sqrt(total(Lp[*,6:7]^2,2,/double))
Fn=sqrt(total(Ln[*,6:7]^2,2,/double))


plot,r,L[*,3],/psym,/xlog,/ylog
oplot,r,Lp[*,3],color=fsc_color('red'),/psym
oplot,r,Ln[*,3],color=fsc_color('blue'),/psym


end
