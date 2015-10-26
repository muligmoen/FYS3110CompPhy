
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import subprocess


def get_output():
  command = ['./project4']
  res = subprocess.check_output(command,universal_newlines=True)
  return res

simulation = get_output().split('\n')

Lx = int(simulation[1].split()[0])
Ly = int(simulation[2].split()[0])
beta = float(simulation[3].split()[0])


def get_state(N):
  """N is the first line of the state input
  """
  state = np.empty((Ly, Lx), dtype=np.int8)
  for i in range(Ly):
    string = simulation[N+i].split()
    for j in range(Lx):
      state[i,j] = int(string[j])
  return state

istate = get_state(5) # initial state

#This updates the state
def update(data):
  mat.set_data(data)
  return mat

#This feeds update()
def data_gen():
  Nstart = 5
  Nstep = Ly + 1
  N = Nstart + Nstep # The initial is already plotted
  while N + Nstep < len(simulation):
    yield get_state(N)
    N += Nstep


fig, ax = plt.subplots()
mat = ax.matshow(istate, cmap='coolwarm')
ax.axis('off')
ax.set_title(r'Ising 2D-lattice, $\beta = {}, L = {}\times {}$'.format(beta, Lx,Ly))

ani = animation.FuncAnimation(fig, update, data_gen, interval=500,
                              save_count=50,repeat=False)
plt.show()
