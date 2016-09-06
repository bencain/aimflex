#!/usr/bin/env python
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

delta = 0.025
x = y = np.arange(-3.0, 3.0, delta)
X, Y = np.meshgrid(x, y)
Z1 = mlab.bivariate_normal(X, Y, 1.0, 1.0, 0.0, 0.0)
Z2 = mlab.bivariate_normal(X, Y, 1.5, 0.5, 1, 1)
Z = Z2 - Z1  # difference of Gaussians

plt.ion()

f, (ax1, ax2, ax3) = plt.subplots(1,3)

ax1.imshow(Z1, interpolation='bilinear', cmap=cm.RdYlGn,
           origin='lower', extent=[-3, 3, -3, 3],
           vmax=abs(Z).max(), vmin=-abs(Z).max())

ax2.imshow(Z2, interpolation='bilinear', cmap=cm.RdYlGn,
           origin='lower', extent=[-3, 3, -3, 3],
           vmax=abs(Z).max(), vmin=-abs(Z).max())

ax3.imshow(Z, interpolation='bilinear', cmap=cm.RdYlGn,
           origin='lower', extent=[-3, 3, -3, 3],
           vmax=abs(Z).max(), vmin=-abs(Z).max())

[a.set_xticks([]) for a in [ax1,ax2,ax3] ]
[a.set_yticks([]) for a in [ax1,ax2,ax3] ]

ax1.set_title('Data')
ax2.set_title('Fit')
ax3.set_title('Residual')

f.savefig('tmp.png',bbox_inches='tight', transparent=True)


