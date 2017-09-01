import numpy as np
import matplotlib.pyplot as plt


N = np.genfromtxt('project1b_n_1000.txt',skip_header=2) 	 # read file, ignore first two lines
x = N[0:,0]						 							 # x from file
u = N[0:,1] 					 							 # analytical result, extract from file
v = N[0:,2]													 # numerical result, extract from file

plt.plot(x,u)
plt.plot(x,v)
plt.xlabel('x')
plt.ylabel('u(x)')
plt.legend(['analytical','numerical'])
plt.show()