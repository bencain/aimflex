pro test_newaim

ngrid=501L
img=dblarr(ngrid,ngrid)

t=dcomplex(dindgen(ngrid^2) mod ngrid,$
           lindgen(ngrid^2)/ngrid) - dcomplex(ngrid/2d,ngrid/2d)

s=1d
a=20d
n=0.5d

te=1000d
r=1.1d*te
pa=0d
kappa= 0.5d*te/r
gamma=-0.5d*te/r
fflex=-0.5d*te/r^2
gflex= 1.5d*te/r^2

shear=[cos(2d*pa),sin(2d*pa)]*gamma/(1d - kappa)
psi_1=[cos(1d*pa),sin(1d*pa)]*fflex/(1d - kappa)
psi_3=[cos(3d*pa),sin(3d*pa)]*gflex/(1d - kappa)

print,shear
print,psi_1
print,psi_3


bkg=100d

p=[s,a,n,shear,psi_1,psi_3]

I=new_aim(t,p,bkg=bkg)
img[*]=I[*]

img-=min(img)
img/=max(img)

tv,img*255d


end
