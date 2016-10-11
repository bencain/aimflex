import aimflex
import numpy as np
# from matplotlib import pyplot
# from astropy.io import fits
# from astropy.modeling import fitting
# from astropy.io import ascii

im = '../../flexion_development/aimdata/hlsp_frontier_hst_acs-30mas-selfcal_abell2744_f814w_v1.0-epoch2_drz.fits'
wt = '../../flexion_development/aimdata/hlsp_frontier_hst_acs-30mas-selfcal_abell2744_f814w_v1.0-epoch2_wht.fits'
cat= '../../flexion_development/aimdata/test.cat'
# cat= 'aimdata/test.cat'


m = aimflex.AIM()

dims = (11,11)
x,y =np.meshgrid(np.linspace(-(dims[0]-1)/2,(dims[0]-1)/2,dims[0]),
					 np.linspace(-(dims[1]-1)/2,(dims[1]-1)/2,dims[1]))
p = np.zeros((3,3))
p[1,1]=1

data = m(x,y,p)

# print aimflex.finite_diff_hessian(m,data,x,y,p,step=1e-6)

import astropy.modeling as am
m = am.functional_models.Sine1D()

A=2.
w=4.
p=1.2

m.amplitude=A
m.frequency=w
m.phase=p


x=np.arange(0,2*np.pi,0.01)

data = m(x)

print aimflex.finite_diff_hessian(m,data,x,step=1e-9)
print aimflex.finite_diff_hessian(m,data,x,step=1e-8)
print aimflex.finite_diff_hessian(m,data,x,step=1e-7)
print aimflex.finite_diff_hessian(m,data,x,step=1e-6)

h=np.zeros((3,3))

h[0,0] = np.sum(2*np.sin(w*x+p)**2)
h[1,0] = np.sum(2*np.sin(w*x+p)*np.cos(w*x+p) + 2*A*np.sin(w*x+p)*w*np.cos(w*x+p) - 2*data*w*np.cos(w*x+p))
h[0,1] = h[1,0]
h[2,0] = np.sum(4*A*np.sin(w*x+p)*np.cos(w*x+p) - 2*data*np.cos(w*x+p))
h[0,2] = h[2,0]
h[1,1] = np.sum(2*A**2*w**2*(np.cos(w*x+p)**2 - np.sin(w*x+p)**2))
h[1,2] = np.sum(2*A**2*w*np.cos(w*x+p)**2 - 2*A**2*w*np.sin(w*x+p)**2 + 2*data*A*w*np.sin(w*x+p))
h[2,1] = h[1,2]
h[2,2] = np.sum(2*A**2*np.cos(w*x+p)**2 - 2*A**2*np.sin(w*x+p)**2 + 2*data*A*np.sin(w*x+p))


print h
# aimflex.fit_dataset(im, wt, cat, 'test_output.txt', rscale=3.,
# 					ntag='NUMBER',xtag='X_IMAGE',ytag='Y_IMAGE',atag='A_IMAGE',
# 					outdir='../../flexion_development/out/',save_fig=True)

