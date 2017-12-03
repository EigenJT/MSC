import numpy as n
import pyspec as p
import time as ti
import random as r
import numpy.linalg as lin
from scipy.optimize import minimize, curve_fit
from scipy.special import erf
import matplotlib.pyplot as plt
from sympy.physics.wigner import wigner_6j

prefilename = 'run05'
postfilename = '_plot.txt'
filenums = ['094','095','096','097','098','099',100,101,102,103,107,108,109,110,111,112,113,114,115,116,117,118,119,120]

for y in n.arange(0,len(filenums),1):
	filenums[y] = str(filenums[y])

n_files = 0
data_list = []
for x in filenums:
	num = str(x)
	filename=prefilename+num+postfilename
	f = open(filename,'r')
	data = n.loadtxt(f)
	data_list.append(data)
	n_files = n_files+1


file_combs = [['094','095'],['096','097'],['098','099'],[100,101],[102,103],[107,108],[109,110],[111,112,],[113,114],[115,116],[117,118],[119,120],[121,122]]
powers = [24,18.9,8.7,4.5,0.79,22.1,28.1,36.4,39.1,50.1,108,12.5,24.4]
nms = []

for y in n.arange(0,len(file_combs),1):
	full_data = n.append(data_list[y],data_list[y+1],axis = 0)
	full_data[:,1] = full_data[n.argsort(full_data[:,0]),1]
	full_data[:,0] =  full_data[n.argsort(full_data[:,0]),0]
	name = str(file_combs[y][0])+'_'+str(file_combs[y][1])
	nms.append(name)
	f = open(name+'.txt','w')
	for x in full_data:
		f.write('%f  %f\n'%(float(x[0]),float(x[1])))
	f.close()

nms_index = 0
for power in powers:
	plt.clf()
	power = power*1e-3
	vec[15] = nms[nms_index]+'.txt'
	vec[12] = power
	execfile('OPV4.py')
	nms_index = nms_index+1
	fig_name = nms[nms_index]+'.png'
	plt.savefig(fig_name)
	
