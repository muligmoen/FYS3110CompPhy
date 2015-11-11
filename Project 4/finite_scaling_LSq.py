import numpy as np
import matplotlib.pyplot as plt

x=np.array([20,40,60,80])		# L values
y=np.array([-1,0.2,0.9,2.1])	# corresponding approximated T_c values

A = np.vstack([x, np.ones(len(x))]).T 	# setting up matrix system for least squares method
a, c = np.linalg.lstsq(A, y)[0]			# solving the system, obtaining the slope

plt.plot(x, y, 'o', label='$T_c (L)$', markersize=10)
plt.plot(x, a*x + c, 'r', label='Best fit')
plt.title('Fitting line through $T_c(L)$ to obtain the slope $a=$ {}'.format(a))
plt.xlabel('L')
plt.ylabel('$T_c (L)$')
plt.legend()
plt.show()