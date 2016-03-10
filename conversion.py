from sympy import *

E,g,F,G,S,T,U,V = symbols("E g F G S T U V")

source(solve)

# print solveset([Eq(S,(E+g)/(1+E*conjugate(g))),
# 		 		Eq(T,(F-conjugate(E)*G)/(1+conjugate(E)*g)),
# 		 		Eq(U,(F-E*conjugate(F))/(1+E*conjugate(g))),
# 		 		Eq(V,(G-E*F)/(1+E*conjugate(g))-V)],
# 		 		[E,g,F,G])

print solve( Eq(37,3*g - 4),g)

# import numpy as np
# import matplotlib.pyplot as plt
# from scipy.stats import gaussian_kde
# 
# def stuv(e1,e2,g1,g2,F1,F2,G1,G2):
# 	E=e1+1j*e2
# 	g=g1+1j*g2
# 	F=F1+1j*F2
# 	G=G1+1j*G2
# 	
# 	M=1. + E*np.conj(g)
# 
# 	S=(E+g)/M
# 	T=(F-np.conj(E)*G)/np.conj(M)
# 	U=(F-E*np.conj(F))/M
# 	V=(G-E*F)/M
# 	
# 	return S.real,S.imag,\
# 		   T.real,T.imag,\
# 		   U.real,U.imag,\
# 		   V.real,V.imag,\
# 
# 
# h=1e-4
# N=1000
# 
# 
# E=2*(np.random.rand(N) + 1j*np.random.rand(N)) - (1. + 1j*1.)
# E=E[abs(E)<1.][0:N]
# 
# g=2*(np.random.rand(N) + 1j*np.random.rand(N)) - (1. + 1j*1.)
# F=0.1*(2*(np.random.rand(N) + 1j*np.random.rand(N)) - (1. + 1j*1.))
# G=0.1*(2*(np.random.rand(N) + 1j*np.random.rand(N)) - (1. + 1j*1.))
# 
# 
# mats=[]

	




# E*=0.7
# g*=2.
# F*=0.1
# G*=0.1
# 
# 
# S,T,U,V
# 
# 
# Sxy=np.vstack([S.real,S.imag])
# Txy=np.vstack([T.real,T.imag])
# Uxy=np.vstack([U.real,U.imag])
# Vxy=np.vstack([V.real,V.imag])
# 
# Sz = gaussian_kde(Sxy)(Sxy)
# Tz = gaussian_kde(Txy)(Txy)
# Uz = gaussian_kde(Uxy)(Uxy)
# Vz = gaussian_kde(Vxy)(Vxy)
# 
# fig, ((axS,axT),(axU,axV)) = plt.subplots(2,2)
# axS.scatter(S.real, S.imag, c=Sz, s=100, edgecolor='')
# axT.scatter(T.real, T.imag, c=Tz, s=100, edgecolor='')
# axU.scatter(U.real, U.imag, c=Uz, s=100, edgecolor='')
# axV.scatter(V.real, V.imag, c=Vz, s=100, edgecolor='')
# 
# plt.show()