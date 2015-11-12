import numpy as np
import matplotlib.pyplot as plt

# Fitting the constant a for finite scaling using least squares method.

L = np.array([20,40,60,80])			# L values
T_C = np.array([2.34,2.31,2.32,2.27])	# corresponding approximated T_c(L) values

a = np.zeros(6)

a[0] = (T_C[3]-T_C[2])/(L[3]-L[2])
a[1] = (T_C[3]-T_C[1])/(L[3]-L[1])
a[2] = (T_C[3]-T_C[0])/(L[3]-L[0])
a[3] = (T_C[2]-T_C[1])/(L[2]-L[1])
a[4] = (T_C[2]-T_C[0])/(L[2]-L[0])
a[5] = (T_C[1]-T_C[0])/(L[1]-L[0])

print a

a = sum(a)/6

print 'Mean value of a: ', a

print 'Finite scaled approximation of T_c: ', T_C[3]-a*L[3]**-1