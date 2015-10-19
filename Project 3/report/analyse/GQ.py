import matplotlib.pyplot as plt

exact = 0.192765711

# limit = 3
def rel_error(x):
  return abs(x-exact)/exact

N = [10, 15, 20, 25, 30, 35]

l3 = [0.071979677, 0.239088291, 0.156139391, 0.195816565,
      0.177282956, 0.189923139]
l5 = [0.011164658, 0.315862997, 0.096788762, 0.240135013,
      0.146370604, 0.204704919]
l10 = [0.000011218, 0.135901531, 0.009799632, 0.307075670, 
       0.051521493, 0.298264908]

laguerre = [0.186457345, 0.189758982, 0.191081780,
            0.191740740, 0.192113712, 0.192343301]


n = len(N);
rel_err = [[0]*n, [0]*n, [0]*n, [0]*n]
for ii in range(n):
  rel_err[0][ii] = rel_error(l3[ii])
  rel_err[1][ii] = rel_error(l5[ii])
  rel_err[2][ii] = rel_error(l10[ii])
  rel_err[3][ii] = rel_error(laguerre[ii])
  
plt.plot(N, rel_err[0], linewidth=2, label='Legendre, limit = 3')
plt.plot(N, rel_err[1], linewidth=2, label='Legendre, limit = 5')
plt.plot(N, rel_err[2], linewidth=2, label='Legendre, limit = 10')
plt.plot(N, rel_err[3], linewidth=2, label='Laguerre')

plt.xlabel('N')
plt.ylabel('Relative error')
plt.xlim((N[0], N[-1]))
plt.ylim((0, 1))
plt.legend()

plt.savefig('GQ_error.png',dpi=400, bbox_inches='tight')

plt.show()

