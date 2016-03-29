import aim
import numpy as np
from matplotlib import pyplot

# Randomly select a set of variables

E = (2*np.random.random(100) - 1.) + 1j*(2*np.random.random(100) - 1.)

elim=0.7
E = E[abs(E)<elim][0] # pick the first one in the allowable range

g = (2*np.random.random() - 1.) + 1j*(2*np.random.random() - 1.)
G = (2*np.random.random() - 1.) + 1j*(2*np.random.random() - 1.)
F = (2*np.random.random() - 1.) + 1j*(2*np.random.random() - 1.)

flexmax=0.1

G *= flexmax
F *= flexmax

print E,g,G,F

S,T,U,V = aim.STUV(E,g,G,F)

print S,T,U,V

Etest = (2*np.random.random(1000) - 1.) + 1j*(2*np.random.random(1000) - 1.)
Etest = Etest[abs(Etest)<0.99][0:100] # pick the first one in the allowable range

gtest = (S*(1.+Etest*np.conj(Etest)) - Etest*(1.+S*np.conj(S)))/(1.+S*Etest*np.conj(S*Etest))

Ftest = ((T+S*np.conj(V))/(1.-S*np.conj(S))) - np.conj(gtest)*((V+S*T)/(1.-S*np.conj(S)))
Gtest = ((V+S*T)/(1.-S*np.conj(S))) - gtest*((T+np.conj(S)*V)/(1.-S*np.conj(S)))

# print abs(E - Etest)
# print abs(g - gtest)
# print abs(F - Ftest)
# print abs(G - Gtest)

a = (T*U-V*np.conj(U))/(1.-S*np.conj(S))
b = (Ftest**2 - np.conj(Ftest)*Gtest)/(1.-gtest*np.conj(gtest))

pyplot.plot(abs(E - Etest),abs(g - gtest))
pyplot.plot(abs(E - Etest),abs(F - Ftest))
pyplot.plot(abs(E - Etest),abs(G - Gtest))
pyplot.plot(abs(E - Etest),abs(b - a))

pyplot.show()
