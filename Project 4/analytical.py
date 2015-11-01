import matplotlib.pyplot as plt
import numpy as np



"""
Numerical solution of the 2x2 Ising system using numerical approximations
"""

def Z(beta):
  return 12 + 4*np.cosh(8*beta)

def lnZ(beta):
  return np.log(Z(beta))

T = np.linspace(1,5,30)
beta = T**-1
Z_beta = Z(beta)

M = (4*2*np.exp(8*beta) + 2*8)/Z_beta/4
E = (-8*2*np.exp(8*beta) + 8*2*np.exp(-8*beta))/Z_beta/4

M2 = ((4**2)*2*np.exp(8*beta) + (2**2)*8)/Z_beta/(4*4)
E2 = ((-8)**2*2*np.exp(8*beta) + (8**2)*2*np.exp(-8*beta))/Z_beta/(4*4)

"""
plt.plot(T, M,label='Magnetisation',linewidth=2)
plt.plot(T, M2, label='Magnetisation squared',linewidth=2)

plt.xlabel('Temperature T')

plt.legend()
plt.show()



plt.plot(T, E, label='Energy',linewidth=2)
plt.plot(T, E2, label='Energy squared', linewidth=2)


plt.xlabel('Temperature T')

plt.legend()
plt.show()
"""

Cv = beta*beta/4*(E2-E*E)
chi = beta*4*(M2-M*M)

plt.plot(T, Cv, label=r'$C_v$',linewidth=2)
plt.plot(T, chi, label=r'$\chi$',linewidth=2)

plt.xlabel('Temperature T')

plt.legend()
plt.show()