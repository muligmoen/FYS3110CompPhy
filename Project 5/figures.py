import subprocess
import matplotlib.pyplot as plt
import numpy as np

def get_stdout(cmd):
  p = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
  return p.stdout.read()

def plot_stuff():
  fig, ax = plt.subplots()

  N = len(u_ANA)
  x = np.linspace(0, 1, N)
  plt.plot(x, u_ANA, linewidth=2, label='Analytical solution')
  plt.plot(x, u_eulerF, linewidth=2, label='Euler forward')
  plt.plot(x, u_eulerB, linewidth=2, label='Euler backward')
  plt.plot(x, u_CN, linewidth=2, label='Crank Nicolson')
  plt.plot(x, u_MC, linewidth=2, label='Monte Carlo')
  plt.plot(x, u_MCG, linewidth=2, label='Monte Carlo variable timestep')

  plt.title('Diffusion of particles',fontsize=16)
  plt.xlabel('x')
  plt.ylabel('Concentration')
  ax.axvspan(-0.2, 0, ymin=0, ymax=1, alpha=0.5, color='red')
  ax.axvspan(1, 1.2, ymin=0, ymax=1, alpha=0.5, color='red')

  ax.text(-0.1, 0.6, r'axon (presynaptic)', rotation='vertical', fontsize=14)
  ax.text(1.1, 0.65, r'dendrite (postsynaptic)', rotation='vertical', fontsize=14)

  plt.axis([-0.2, 1.2, 0, 1])

  plt.legend(loc='upper right')
  plt.show()
  
  
alpha = 1.0/2 # Stability criteria
dx = 1/100
dt = dx*dx*alpha
T = 0.1

cmd = ['./project5', str(dt), str(dx), str(T)]

output = get_stdout(cmd)

lines = output.split('\n')

u_ANA = list(map(float, lines[0].strip('[]').split()))
u_eulerF = list(map(float, lines[1].strip('[]').split()))
u_eulerB = list(map(float, lines[3].strip('[]').split()))
u_CN = list(map(float, lines[5].strip('[]').split()))
u_MC = list(map(float, lines[7].strip('[]').split()))
u_MCG = list(map(float, lines[10].strip('[]').split()))

E = [float(lines[2]), float(lines[4]), float(lines[6]), float(lines[8]), float(lines[11])]
print(E)

plot_stuff()