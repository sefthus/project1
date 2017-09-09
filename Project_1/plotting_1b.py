import numpy as np
import matplotlib.pyplot as plt


N1 = np.genfromtxt('project1b_n_10.txt',skip_header=2) 	 # read file, ignore first two lines
x1 = N1[0:,0]						 # x from file
u1 = N1[0:,1] 					 # analytical result, extract from file
v1 = N1[0:,2]				 # numerical result, extract from file

N2 = np.genfromtxt('project1b_n_100.txt',skip_header=2) 	 # read file, ignore first two lines
x2 = N2[0:,0]						 # x from file
u2 = N2[0:,1] 					 # analytical result, extract from file
v2 = N2[0:,2]				 # numerical result, extract from file

N3 = np.genfromtxt('project1b_n_1000.txt',skip_header=2) 	 # read file, ignore first two lines
x3 = N3[0:,0]						 # x from file
u3 = N3[0:,1] 					 # analytical result, extract from file
v3 = N3[0:,2]				 # numerical result, extract from file

def plot_formatting():
	"""pretty formatting of plots """
	plt.rc('text',usetex=True) # only works with matplotlib version below 1.5.2, comment out elsewise
	axis_font={'family': 'serif','serif':'Computer Modern Roman','size':16}
	plt.rc('font',**axis_font)
	plt.rc('font',weight ='bold')
	plt.xticks(fontsize=16)
	plt.yticks(fontsize=16)

plot_formatting()
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')

plt.subplot(3,1,1)
plt.plot(x1,u1,x1,v1,linewidth=2.0)
plt.legend(['Analytical','Numerical'],loc ='upper right')
plt.setp(plt.gca().get_legend().get_texts(), fontsize=16)
plt.title('n=10',fontsize=16)
plt.xlabel('$x$')
plt.ylabel('$u(x)$')

plt.subplot(3,1,2)
plt.plot(x2,u2,x2,v2,linewidth=2.0)
plt.legend(['Analytical','Numerical'],loc ='upper right')
plt.setp(plt.gca().get_legend().get_texts(), fontsize=16)
plt.title('n=100',fontsize=16)
plt.xlabel('$x$')
plt.ylabel('$u(x)$')

plt.subplot(3,1,3)
plt.plot(x3,u3,x3,v3,linewidth=2.0)
plt.legend(['Analytical','Numerical'],loc ='upper right')
plt.setp(plt.gca().get_legend().get_texts(), fontsize=16)
plt.title('n=1000',fontsize=16)
plt.xlabel('$x$')
plt.ylabel('$u(x)$')


plt.suptitle('Exact and numerical solution of set of linear equations',fontsize=16)
plt.tight_layout()
plt.subplots_adjust(top=0.88)
plt.show()
	
