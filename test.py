import aimflex
import numpy as np
import corner

from matplotlib import pyplot
from matplotlib.ticker import MaxNLocator

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
m.F1_0 = -0.005
m.G1_0 = -0.001

dims = (51,51)
x,y =np.meshgrid(np.linspace(-(dims[0]-1)/2,(dims[0]-1)/2,dims[0]),
					 np.linspace(-(dims[1]-1)/2,(dims[1]-1)/2,dims[1]))
p = np.zeros((3,3))
p[1,1]=1

true_params = np.copy(m.parameters)

true_img = aimflex.window_image(m(x,y,p))
data = aimflex.window_image(m(x,y,p) + (2*np.random.random(dims) - 1.))
weights = aimflex.window_image(np.ones(dims))


aimflex.set_gaussian_pars(data,m,weights=weights)
aimflex.set_limits(data,m)

samples = aimflex.fit_image(m,data,weights,p,verbose=True)

m.parameters = samples[-1,:]

last = aimflex.window_image(m(x,y,p))

fig = corner.corner(samples, labels=m.param_names,
					truths=true_params)
fig.savefig("aim-triangle.png")


f, axarr = pyplot.subplots(2,2)
axarr[0,0].imshow(data)
axarr[0,1].imshow(true_img)
axarr[1,0].imshow(last)
axarr[1,1].imshow(data-last)

f.savefig('4img.png')

pyplot.clf()
fig, axes = pyplot.subplots(m.parameters.size, 1, sharex=True, figsize=(8, 20))

for i in range(m.parameters.size):
	axes[i].plot(samples[:, i].T, color="k", alpha=0.4)
	axes[i].yaxis.set_major_locator(MaxNLocator(5))
	axes[i].axhline(true_params[i], color="b", lw=2)
	axes[i].set_ylabel(m.param_names[i])

fig.tight_layout(h_pad=0.0)
fig.savefig("aim-time.png")
