import subprocess
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

def get_stdout(cmd):
  p = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
  return p.stdout.read()

command = ['./project5']

strings = get_stdout(command).split('\n')


N = int(strings.pop(0))
Nt = int(strings.pop(0))

u = [[0]*N]*Nt

for ii in range(Nt):
  arr = strings[ii].split()[1:-1]
  u[ii] = [float(x) for x in arr]

fig, ax = plt.subplots()

x = np.linspace(0, 1, N)
line, = ax.plot(x, u[0], linewidth=2)
plt.xlabel('x')
plt.ylabel('Concentration')

plt.ylim(-1, 1)

#plt.ylim(0, 150) #MC with 150 particles ?


def animate(i):
  line.set_ydata(u[i])
  return line,



ani = animation.FuncAnimation(fig, animate, range(Nt), interval=250)

plt.show()
#for i in range()
#print(strings)