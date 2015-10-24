
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import subprocess


def runProcess(exe):    
    p = subprocess.Popen(exe, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    while(True):
      retcode = p.poll() #returns None while subprocess is running
      line = p.stdout.readline()
      yield line
      if(retcode is not None):
        break


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
    #N += 1
  return state

istate = get_state(5)

#plt.matshow(istate)
#plt.show()
  

def update(data):
  mat.set_data(data)
  return mat

def data_gen():
  Nstart = 5
  Nstep = Ly + 1
  N = Nstart + Nstep
  while N + Nstep < len(simulation):
    yield get_state(N)
    N += Nstep


fig, ax = plt.subplots()
mat = ax.matshow(istate, cmap='coolwarm')
ax.axis('off')

ani = animation.FuncAnimation(fig, update, data_gen, interval=500,
                              save_count=50,repeat=False)
plt.show()
"""

def line_numbers():
  N = 5 + Lx + 2
  while N < len(simulation):
    N += Lx + 2
    yield N
    
def yield_state():
  yield get_state(line_numbers)
  
def init():
  mat = ax.matshow(istate)
  return mat

def update(data):
  mat.set_data(data)
  return mat
  
fig, ax = plt.subplots()
mat = ax.matshow(istate)
  
ani = animation.FuncAnimation(fig, update, yield_state,
                              interval=500, save_count=50)
plt.show()
"""
