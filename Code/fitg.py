import numpy as n
import pyspec as p
import time as ti
import random as r
import numpy.linalg as lin
from scipy.optimize import minimize, curve_fit
from scipy.special import erf
import matplotlib.pyplot as plt
from sympy.physics.wigner import wigner_6j




#------------------------INTRO----------------------------

print '\
		                                   _ _ _ _ _ _ _            \n\
		         _ _ _ _ _       /        /                         \n\
		        /               /        /                          \n\
		       /               /        /                           \n\
		      /_ _ _ _   - - -/- - -   /                            \n\
		     /         o     /        /                             \n\
		    /         /     /        /        _ _ _ _               \n\
		   /         /     /        /        /      /               \n\
		  /         /     /        /               /                \n\
		                          /_ _ _ _ _ _ _ _/                 \n\
		                                                            \n\
		                                                            \n\
		                                                 JRT(2016)  \n\
		                                                            '

print "\n\nThis is a general purpose fitting routine for nuclear spectra.\n\n\
Here are a few instructions:\n\n\
1) You'll be asked the name of the data file. The data should be a .txt\n\
file with the first column as the frequency and the second column as the counts.\n\
It should also be located in the same directory as this script.\n\n\
2) Next, you'll be asked for the nuclear spin (I) of whatever isotope you're fitting.\n\n\
3) You will then be asked for the spins of the upper J state and lower J state.\n\n\
4) Once you've input all the information, the data will be plotted. \n\n\
5) Now, some initial guesses would be useful to help the fitting. \n\
Bradly Cheal's 'viewer.py' program is a god send for this step. Seriously, use it.\n\n\
6) Once you've entered those guesses, the best fit will be found and plotted. \n\n\
7) This program will have a few outputs. It will have a fit spectrum, with the\n\
hyperfine splitting parameters, centroid and reduced chi squared on it.\n\
Additionally, an output file will be generated with the fit parameters written to them if\n\
you need to replot the spectrum or store the results.\n"


#--------------Information input-------------------------

c = raw_input('Proceed?(y/n)')#getting the information about the isotope
if c in ['n','N','No','n']:
	print 'Awwww...'
	exit()
elif c in ['y','Y','Yes','yes']:
	print 'Please choose the data file.'
	filename = raw_input('Filename:')
	data = n.loadtxt(open(filename,"r"))
	freq = data[:,0]#import frequency
	cts = data[:,1]#import counts
	err = n.sqrt(cts)
	print 'Correcting for zero error...'
	for ch in n.arange(0,len(cts),1):#if anything is zero, then the error is set to one. 
		if err[ch] == 0:
			print 'Found a zero!'
			err[ch] = 1.0
	print 'Done'
	I = float(raw_input('Nuclear Spin (I):'))
	J_u = float(raw_input('Upper J state (J_u):'))
	J_l = float(raw_input('Lower J state (J_l):'))
	plt.ion()
	plt.errorbar(freq,cts,err,fmt='k.',label = filename)
	plt.xlabel('Frequency (MHz)')
	plt.ylabel('Counts')
	plt.legend()
	plt.draw()



#------------------------Definitions--------------------

def F(I,J_u,J_l):#returns all the possible transistions Fu->Fl J_u, J_l
	J_us = n.arange(-J_u,J_u+1,1)
	J_ls = n.arange(-J_l,J_l+1,1)
	F_us = I + J_us
	F_ls = I + J_ls
	trans = []
	ntrans = 0
	for u in F_us:
		for l in F_ls:
			ch = n.abs(u-l)
			if ch == 1:
				trans.append(u)
				trans.append(l)
				trans.append(J_u)
				trans.append(J_l)
				ntrans = ntrans + 1
			if ch == 0:
				if u != 0:
					trans.append(u)
					trans.append(l)
					trans.append(J_u)
					trans.append(J_l)
					ntrans = ntrans + 1
	Ftemp = n.array(trans)
	F = Ftemp.reshape((ntrans,4))
	return F

def voigt(xdata,amp,cent,sig,ep):#define the voigt peak. If ep is less than 0 or greater than 1, mess up the profile
    x = xdata
    C = cent
    S = sig
    A = amp
    E = ep
    Sg = S/n.sqrt(2*n.log(2))
    vmodel = A*(1.0-E)/(Sg*n.sqrt(2.0*n.pi))*n.exp(-(x-C)**2/(2.0*Sg**2))+ (E*A)/(n.pi)*(S)/((x-C)**2+S**2)
    if 0>E or 1<E:
    	vmodel = vmodel + 50000
    return vmodel



def chi2(x0):
	nu = cts.size - len(x0) - 1 #degrees of freedom
	if x0[2] == 0:
		nu = nu + 1
	if x0[3] == 0:
		nu = nu + 1
	model = spec(freq,x0[0],x0[1],x0[2],x0[3],x0[4],x0[5],x0[6],x0[7],x0[8])
	chi2 = (nu**(-1))*n.sum(((cts-model)*(err**(-1)))**(2))
	return chi2

def chi2r(x0):
	nu = cts.size - len(x0) - 1 #degrees of freedom
	if x0[1] == 0:
		nu = nu + 1
	if x0[2] == 0:
		nu = nu + 1
	model = specr(freq,x0[0],x0[1],x0[2],x0[3],x0[4],x0[5],x0[6],x0[7])
	chi2 = (nu**(-1))*n.sum(((cts-model)*(err**(-1)))**(2))
	return chi2


#Define K
def K(F,I,J):
	Kk = F*(F+1.0) - I*(I+1.0) - J*(J+1.0)
	return Kk

#Define Beta
def Beta(K,I,J):
	if I <= 0.5 or J <= 0.5:#check for B coeffs. Set to zero if they're not important
		Beta = 0
	else :
		Beta = (3.0*K*(K+1.0)-4.0*I*(I+1.0)*J*(J+1.0))/(8.0*I*(2.0*I-1.0)*J*(2.0*J-1.0))
	return Beta

F = F(I,J_u,J_l)
amps = n.zeros((len(F[:,0]),5))
amps[:,1:5] = F

def pkpos(A_u,A_l,B_u,B_l,cent):#gets the positions of the peaks
	ntrans = len(amps[:,0])
	pos = n.zeros((ntrans,1))
	for p in n.arange(0,ntrans,1):
		F_u = amps[p,1]
		F_l = amps[p,2]
		J_u = amps[p,3]
		J_l = amps[p,4]
		K_u = K(F_u,I,J_u)
		K_l = K(F_l,I,J_l)
		Beta_u = Beta(K_u,I,J_u)
		Beta_l = Beta(K_l,I,J_l)
		pos[p,0] = A_u*K_u/2 + Beta_u*B_u - A_l*K_l/2 - Beta_l*B_l + cent
	peaks = n.zeros((ntrans,6))
	peaks[:,0] = pos[:,0]
	peaks[:,1:6] = amps
	return peaks

def spec(xdat,A_u,A_l,B_u,B_l,cent,Amp,bg,sig,ep):
	num_p = len(amps[:,0])
	pks = pkpos(A_u,A_l,B_u,B_l,cent)
	spc = n.zeros(xdat.shape)
	for pk in n.arange(0,num_p,1):
		spc_tmp = voigt(xdat,Amp*pks[pk,1],pks[pk,0],sig,ep)
		spc = spc + spc_tmp
	spc = spc+bg
	return spc



def lbd(x):
	y = x - 0.01*n.abs(x) - 10**(-20)
	return y
def ubd(x):
	y = x + 0.01*n.abs(x) + 10**(-20)
	return y

def errrat(A_u,A_u_er,A_l,A_l_er):
	err = n.sqrt((A_u_er/A_l)**2+(A_u*A_l_er/(A_l**2))**2)
	return err



#--------------Amplitudes----------------------------#
print 'Assigning transitions to peaks...'
a = 0
for jj in n.arange(0,len(F[:,0]),1):
		amps[jj,0] = (2.0*F[jj,0]+1.0)*(2.0*F[jj,1]+1.0)*(wigner_6j(F[jj,1],F[jj,0],1,F[jj,2],F[jj,3],I))**2

amps[:,0] = amps[:,0]/n.sum(amps[:,0])#normalize the amplitudes to get the ratios
ord = n.argsort(amps[:,0])

amps_temp = n.zeros(amps.shape)
i = 0
for o in ord:
	amps_temp[i,:] = amps[o,:]
	i = i+1
amps = amps_temp#now have the peaks ordered by their theoretical amplitudes
for oo in n.arange(0,len(F[:,0]),1):
	F1 = amps[oo,1]
	F2 = amps[oo,2]
	amp_p = amps[oo,0]
	ti.sleep(0.3)
	print 'F = %f -> F = %f: %f \n'%(n.around(F1,1),n.around(F2,1),amp_p)

print 'Done'

#----------------Guesses-----------------------------#



A_ug = float(raw_input('A_u guess (MHz):'))
A_lg = float(raw_input('A_l guess (MHz):'))
if J_u > 0.5 and I > 0.5:
	B_ug = float(raw_input('B_u guess (MHz):'))
else:
	B_ug = 0
if J_l > 0.5 and I > 0.5:
	B_lg = float(raw_input('B_l guess (MHz):'))
else:
	B_lg = 0
cent_g = float(raw_input('Centroid guess (MHz):'))


#------------Fitting----------------#
rat = raw_input('Would you like to use a ratio of A_u/A_l?(y/n)')
sran = n.arange(n.min(freq),n.max(freq),0.1)

if rat in ['y']:
	rat_in = float(raw_input('A_u/A_l?'))
	
	def specr(xdat,A_u,B_u,B_l,cent,Amp,bg,sig,ep):
		num_p = len(amps[:,0])
		pks = pkpos(A_u,A_u/rat_in,B_u,B_l,cent)
		spc = n.zeros(xdat.shape)
		for pk in n.arange(0,num_p,1):
			spc_tmp = voigt(xdat,Amp*pks[pk,1],pks[pk,0],sig,ep)
			spc = spc + spc_tmp
		spc = spc+bg
		return spc		
	guess = n.array([A_ug,B_ug,B_lg,cent_g,150*max(cts),min(cts),70,0.5])
	popt,pcov = curve_fit(specr,freq,cts,p0=guess,sigma=err,absolute_sigma=True)
	#Once we've gotten the ideal parameters, fit again, but with bounds to get the errors
	bds = ((lbd(popt[0]),lbd(popt[1]),lbd(popt[2]),lbd(popt[3]),lbd(popt[4]),lbd(popt[5]),lbd(popt[6]),lbd(popt[7])),
			(ubd(popt[0]),ubd(popt[1]),ubd(popt[2]),ubd(popt[3]),ubd(popt[4]),ubd(popt[5]),ubd(popt[6]),ubd(popt[7])))
	popt_er,pcov_er = curve_fit(specr,freq,cts,p0=popt,sigma=err,absolute_sigma=True,bounds=bds)
	
	fits = specr(sran,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5],popt[6],popt[7])
	specg = specr(sran,guess[0],guess[1],guess[2],guess[3],guess[4],guess[5],guess[6],guess[7])
	plt.clf()
	sran = n.arange(n.min(freq),n.max(freq),0.1)
	plt.plot(sran,specg,'r-',linewidth = 2.0,label='Initial Guess')
	plt.plot(sran,fits, 'm-',label='Best Fit',linewidth = 2.0)
	plt.errorbar(freq,cts,err,fmt='k.',label = 'Counts')
	plt.xlabel('Frequency (MHz)')
	plt.ylabel('Counts')
	plt.legend()
	plt.show()
	
	print 'Reduced Chi square = %f'%(n.around(chi2r(popt),4))
	print 'A_u = %f +/- %f'%(popt[0], pcov_er[0,0]**(0.5))
	print 'A_l = %f +/- %f'%(popt[0]/rat_in, pcov_er[0,0]**(0.5)/rat_in)
	print 'A_u/A_l = %f (fixed)'%(rat_in)
	print 'B_u = %f +/- %f'%(popt[1], pcov_er[1,1]**(0.5))
	print 'B_l = %f +/- %f'%(popt[2], pcov_er[2,2]**(0.5))
	print 'Centroid = %f +/- %f'%(popt[3], pcov_er[3,3]**(0.5))
	print 'FWHM = %f +/- %f'%(popt[6], pcov_er[6,6]**(0.5))
	print 'ep = %f +/- %f'%(popt[7],pcov_er[7,7])
	print 'Bg = %f +/- %f \n'%(popt[5],pcov_er[5,5])
	
	rs = filename.replace('.txt','_results.txt')
	op = open(rs,'w')
	op.write('Reduced Chi square = %f\n'%(n.around(chi2r(popt),4)))
	op.write('A_u = %f +/- %f\n'%(popt[0], pcov_er[0,0]**(0.5)))
	op.write('A_l = %f +/- %f\n'%(popt[0]/rat_in, pcov_er[0,0]**(0.5)/rat_in))
	op.write('A_u/A_l = %f (fixed)\n'%(rat_in))
	op.write('B_u = %f +/- %f\n'%(popt[1], pcov_er[1,1]**(0.5)))
	op.write('B_l = %f +/- %f\n'%(popt[2], pcov_er[2,2]**(0.5)))
	op.write('Centroid = %f +/- %f\n'%(popt[3], pcov_er[3,3]**(0.5)))
	op.write('FWHM = %f +/- %f\n'%(popt[6], pcov_er[6,6]**(0.5)))
	op.write('ep = %f +/- %f\n'%(popt[7],pcov_er[7,7]))
	op.write('Bg = %f +/- %f \n'%(popt[5],pcov_er[5,5]))
	op.close()
	ft = filename.replace('.txt','.png')
	plt.savefig(ft)
else:
	guess = n.array([A_ug,A_lg,B_ug,B_lg,cent_g,150*max(cts),min(cts),70,0.5])
	popt,pcov = curve_fit(spec,freq,cts,p0=guess,sigma=err,absolute_sigma=True)
	#Once we've gotten the ideal parameters, fit again, but with bounds to get the errors
	bds = ((lbd(popt[0]),lbd(popt[1]),lbd(popt[2]),lbd(popt[3]),lbd(popt[4]),lbd(popt[5]),lbd(popt[6]),lbd(popt[7]),lbd(popt[8])),
			(ubd(popt[0]),ubd(popt[1]),ubd(popt[2]),ubd(popt[3]),ubd(popt[4]),ubd(popt[5]),ubd(popt[6]),ubd(popt[7]),ubd(popt[8])))
	popt_er,pcov_er = curve_fit(spec,freq,cts,p0=popt,sigma=err,absolute_sigma=True,bounds=bds,method='dogbox')
	
	fits = spec(sran,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5],popt[6],popt[7],popt[8])
	specg = spec(sran,guess[0],guess[1],guess[2],guess[3],guess[4],guess[5],guess[6],guess[7],guess[8])
	plt.clf()
	plt.plot(sran,specg,'r-',linewidth = 2.0,label='Initial Guess')
	plt.plot(sran,fits, 'm-',label='Best Fit',linewidth = 2.0)
	plt.errorbar(freq,cts,err,fmt='k.',label = 'Counts')
	plt.xlabel('Frequency (MHz)')
	plt.ylabel('Counts')
	plt.legend()
	plt.show()
	
	print 'Reduced Chi square = %f'%(n.around(chi2(popt),4))
	print 'A_u = %f +/- %f'%(popt[0], pcov_er[0,0]**(0.5))
	print 'A_l = %f +/- %f'%(popt[1], pcov_er[1,1]**(0.5))
	print 'A_u/A_l = %f +/- %f'%(popt[0]/popt[1],errrat(popt[0],pcov_er[0,0]**(0.5),popt[1],pcov_er[1,1]**(0.5)))
	print 'B_u = %f +/- %f'%(popt[2], pcov_er[2,2]**(0.5))
	print 'B_l = %f +/- %f'%(popt[3], pcov_er[3,3]**(0.5))
	print 'Centroid = %f +/- %f'%(popt[4], pcov_er[4,4]**(0.5))
	print 'FWHM = %f +/- %f'%(popt[7], pcov_er[7,7]**(0.5))
	print 'Bg = %f +/- %f \n'%(popt[6],pcov_er[6,6])
	print 'ep = %f +/- %f'%(popt[8],pcov_er[8,8])
	
	rs = filename.replace('.txt','_results.txt')
	op = open(rs,'w')
	op.write('Reduced Chi square = %f\n'%(n.around(chi2(popt),4)))
	op.write('A_u = %f +/- %f\n'%(popt[0], pcov_er[0,0]**(0.5)))
	op.write('A_l = %f +/- %f\n'%(popt[1], pcov_er[1,1]**(0.5)))
	op.write('A_u/A_l = %f +/- %f\n'%(popt[0]/popt[1],errrat(popt[0],pcov_er[0,0]**(0.5),popt[1],pcov_er[1,1]**(0.5))))
	op.write('B_u = %f +/- %f\n'%(popt[2], pcov_er[2,2]**(0.5)))
	op.write('B_l = %f +/- %f\n'%(popt[3], pcov_er[3,3]**(0.5)))
	op.write('Centroid = %f +/- %f\n'%(popt[4], pcov_er[4,4]**(0.5)))
	op.write('Bg = %f +/- %f \n'%(popt[6],pcov_er[6,6]))
	op.write('FWHM = %f +/- %f\n'%(popt[7], pcov_er[7,7]**(0.5)))
	op.write('ep = %f +/- %f\n'%(popt[8],pcov_er[8,8]))
	op.close()
	ft = filename.replace('.txt','.png')
	plt.savefig(ft)
























