from numpy import *
import matplotlib.pylab as plt

f=open('outputLOL.txt')


lines=f.readlines()

N = int(lines[1])
rhos = array(list(lines[2]))
rho_0 = float(rhos[0])
rho_inf = float(rhos[1])

rho = linspace(rho_0,rho_inf,N)

psi_1 = [float(x) for x in lines[10].split()]
psi_2 = [float(x) for x in lines[11].split()]
psi_3 = [float(x) for x in lines[11].split()]

plt.plot(rho,psi_1)
plt.plot(rho,psi_2)
plt.plot(rho,psi_3)



