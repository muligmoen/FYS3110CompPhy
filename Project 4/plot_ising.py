import matplotlib.pyplot as plt

filename = 'test.txt'

with open(filename) as f:
  N = int(f.readline().split()[0]);
  f.readline()
  f.readline()
  f.readline()
  f.readline()
  f.readline()
  f.readline()
  f.readline()
  f.readline()
  f.readline() # Input types
  
  T = [0]*N
  E = [0]*N
  M = [0]*N
  Cv = [0]*N
  chi = [0]*N
  ar = [0]*N
  
  for i in range(N):
    line = f.readline();
    T[i], E[i], M[i], Cv[i], chi[i], ar[i] = map(float, line.split())

f, ax = plt.subplots(2,2,sharex=True)

ax[0, 0].plot(T, E, label='Energy', linewidth=2)
ax[0, 0].set_title('$<E>$')

ax[0, 1].plot(T, M, label='Magnetism', linewidth=2)
ax[0, 1].set_title(r'$<M>$')

ax[1, 0].plot(T, Cv, label=r'$C_v$', linewidth=2)
ax[1, 0].set_title(r'$C_v$')

ax[1, 1].plot(T, chi, label=r'$\chi$', linewidth=2)
ax[1, 1].set_title(r'$\chi$')


#plt.legend()
plt.show()
  
  