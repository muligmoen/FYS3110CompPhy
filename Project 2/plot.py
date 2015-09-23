import numpy as np
import matplotlib.pyplot as plt

f=open('test.txt')


lines=f.readlines()

N = int(lines[0].split()[0])
rhos = lines[1].split()
rho_0 = float(rhos[0])
rho_inf = float(rhos[1])

rho = np.linspace(rho_0,rho_inf,N)


psi_1 = [float(x)**2 for x in lines[9].split()]
psi_2 = [float(x)**2 for x in lines[10].split()]
psi_3 = [float(x)**2 for x in lines[11].split()]

plt.plot(rho,psi_1)
plt.plot(rho,psi_2)
plt.plot(rho,psi_3)

plt.show()


