import aimflex
import numpy as np
from matplotlib import pyplot
# from astropy.io import fits
# from astropy.modeling import fitting
# from astropy.io import ascii

im = '../../flexion_development/aimdata/hlsp_frontier_hst_acs-30mas-selfcal_abell2744_f814w_v1.0-epoch2_drz.fits'
wt = '../../flexion_development/aimdata/hlsp_frontier_hst_acs-30mas-selfcal_abell2744_f814w_v1.0-epoch2_wht.fits'
cat= '../../flexion_development/aimdata/test.cat'
# cat= 'aimdata/test.cat'


m = aimflex.get_AIM()

m.logI = 2.
m.g1_0 = -1./7
m.F1_0 = -0.001
m.G1_0 = -0.001

dims = (51,51)
x,y =np.meshgrid(np.linspace(-(dims[0]-1)/2,(dims[0]-1)/2,dims[0]),
					 np.linspace(-(dims[1]-1)/2,(dims[1]-1)/2,dims[1]))
p = np.zeros((3,3))
p[1,1]=1

data = m(x,y,p) + (2*np.random.random(dims) - 1.)
weights = np.ones(dims)

pyplot.imshow(data)
pyplot.show()

chain = aimflex.fit_image(m,data,weights,p,verbose=True)

for i in range(m.parameters.size):
	print i, len(set(chain[:,:,i].flatten()))
