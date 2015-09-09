#!/usr/bin/env python3
from numpy import *
import subprocess
import matplotlib.pylab as plt

filename = 'test.txt'

N = [10, 100, 1000, 1000]

Alg = ['M','T','S','L']

def analyze_file(filename, N):
	f = open(filename)
	for ii in range(2*N+4):
		f.readline()
	time = float(f.readline().split()[0])					#						  N1  N2  N3  N4	
	err  = float(f.readline().split()[0])					#					M	|				 |
	return time, err										#	matrix format:	T	|				 |
															#					S	|				 |
time = zeros((4,4))											# 					L	|				 |	  
err = zeros((4,4))																	
																				    
for i in range(0,4):
	for j in range(0,4):
		command = ['./project1', str(N[i]), filename, Alg[j]]
		subprocess.call(command)
		time[j,i], err[j,i] = analyze_file(filename, N[i])


log_n = [log(10), log(100), log(1000), log(10000)]

plt.plot(log_n,err[0,:],label="Gaussian elimination")
plt.plot(log_n,err[1,:],label="Thomas algorithm")
plt.plot(log_n,err[3,:],label="LU decomposition")
plt.plot(log_n,err[2,:],label="Sparse matrix decomposition")
plt.xlabel('Log number of mesh points')
plt.ylabel('Log relative error')
plt.title('Max relative error for different algorithms')
plt.legend()
plt.show()

plt.plot(log_n,time[0,:],label="Gaussian elimination")
plt.plot(log_n,time[1,:],label="Thomas algorithm")
plt.plot(log_n,time[3,:],label="LU decomposition")
plt.plot(log_n,time[2,:],label="Sparse matrix decomposition")
plt.xlabel('Log number of mesh points')
plt.ylabel('Effective seconds(Clock cycles/MHz))')
plt.title('Comparing time')
plt.legend()
plt.show()

