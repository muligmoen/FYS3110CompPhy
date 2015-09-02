#!/usr/bin/env python3

import matplotlib.pyplot as plt
import sys
import numpy as np

if (len(sys.argv)>1):
  filename = sys.argv[1]
else:
  print("Please give filename")
  sys.exit(1)
  
f = open(filename);

N = int(f.readline())

b = [0]*N
u = [0]*N
x = np.linspace(0, 1, N)

f.readline()
for i in range(N):
  b[i] = float(f.readline())

f.readline()
for i in range(N):
  u[i] = float(f.readline())

f.close()

plt.plot(x, b, label='numerical')
plt.plot(x, u, label='exact')
plt.legend()
plt.show()
