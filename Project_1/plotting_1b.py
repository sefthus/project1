import numpy as np
import matplotlib.pyplot as plt


N = np.genfromtxt('project1c_n_10.txt',skip_header=2) 	 # read file, ignore first two lines
x = N[0:,0]						 # x from file
u = N[0:,1] 					 # analytical result, extract from file
v = N[0:,2]				 # numerical result, extract from file

def plot_formatting(fam='serif',fam_font='Computer Modern Roman',font_size=16,tick_size=16):
	""" you get to define what font and size of xlabels and axis ticks you"""
	"""like, if you want bold text or not.								  """

	plt.rc('text',usetex=True)
	axis_font={'family': fam,'serif':[fam_font],'size':font_size}
	plt.rc('font',**axis_font)
	plt.rc('font',weight ='bold')
	plt.xticks(fontsize=tick_size)
	plt.yticks(fontsize=tick_size)
#plot_formatting()
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
plt.plot(x,u)
plt.plot(x,v)
plt.xlabel('x')
plt.ylabel('u(x)')
plt.legend(['analytical','numerical'])
plt.show()
	