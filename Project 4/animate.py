import matplotlib.pyplot as plt
import subprocess
import random


def get_output(command):
  res = subprocess.check_output(command,universal_newlines=True)
  return res

def get_magnetisation(N, L, beta):
  command = ['./project4', '-P', 'M', str(N), str(L), str(beta)]
  out = get_output(command).split('\n')[0]
  
  str_M = out.split()
  spins = L*L
  if (float(str_M[-1]) > 0):
    polarity = 1
  else:
    polarity = -1
  pol_spin = polarity/spins
  M = [pol_spin*float(x) for x in str_M]
  return M

def get_energy(N, L, beta):
  command = ['./project4', '-P', 'E', str(N), str(L), str(beta)]
  out = get_output(command).split('\n')[0]
  
  str_E = out.split()
  spins = L*L
  inv_spins = 1.0/spins;
  E = [inv_spins*float(x) for x in str_E]
  return E


N = 20
L = 2
beta = 0.001


plt.subplot(2,1,1)
plt.ylim([-0.2,1.05])
plt.ylabel('Net magnetisation per spin')
plt.xlabel('Number of cycles')

for i in range(10):
  M = get_magnetisation(N, L, beta)
  plt.plot(M,linewidth=2)


plt.subplot(2,1,2)
plt.ylabel('Energy per spin')
plt.xlabel('Number of cycles')
plt.ylim([-4.1,4.1])

for i in range(10):
  E = get_energy(N, L, beta)
  plt.plot(E, linewidth=2)
plt.show()
