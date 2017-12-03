#This version of the code outputs the likelihood of an atom making it to the LCR in the same ground state that it started in.

import numpy as n
import math as m
import cmath as cm
import time as ti
from decimal import *
import random as r
import numpy.linalg as lin
from sympy.physics.wigner import wigner_6j, wigner_3j
from scipy.special import wofz
import matplotlib.pyplot as plt



#----------------------INTRO---------------------------


#------------------Preamble------------------------------------
# Important constants

kB = 8.6173324 * 10**(-5) #eV/K boltzmanns constant
kB_J = 1.38064852e-23#boltzmann constant in J/K
amu = 931.4940954*10**6 #eV/c^2 rest energy of one atomic mass unit
amu_kg = 1.66054e-27 #amu in kg
c = 299792458 #m/s speed of light
alpha = 0.0072973525664 #fine structure constant
m_e = 0.5109989461*10**6#eV/c^2 mass of electron
c_e = 1.6021766208*10**(-19) #charge of electron in Coulombs
r_e = 2.817*10**(-15)#classical electron radius in m
kB_cc = kB/(c**2)#kB/c^2 boltzmann constant per c^2
hbar = 6.582119514*10**(-16) #reduced planck constant in eV*s
h_js = 6.626070040*10**(-34) #planks constant in J s
e_0 = 8.854187817*10**(-6)/(6.242*10**(24)) #vacuum perimittivity in Coulombs^2/eV*m
a_0 =  5.29*10**(-11) #bohr radius in m



def GSP(I,J_l):#Returns a random F state for ground state 
	Js = n.arange(-J_l,J_l + 1,1)
	Fs = Js + I
	F_mg = 2*Fs+1
	tot = n.sum(F_mg)
	F_mg = F_mg/tot
	F_mat = n.zeros([len(Js),2])
	F_mat[:,0] = Fs
	F_mat[:,1] = F_mg
	rnum = r.uniform(0,1)
	lower = 0
	for f in n.arange(0,len(Js)):
		upper = lower + F_mat[f,1]
		if rnum == 1:
			F_state = Fmat[len(Js)-1,0]
		elif lower <= rnum < upper:
			F_state = F_mat[f,0]
		lower = upper
	return F_state



def Boltz(M,T):#returns random z-velocity from a Boltzmann distribution (m/s)
	sig = (kB*T/(2.0*n.pi*M/c**2))**(0.5)
	ran1 = n.random.uniform(0,1)
	ran2 = n.random.uniform(0,1)
	vel = sig*n.sqrt(-2.0*n.log(ran1))*n.cos(2.0*n.pi*ran2)
	return vel

		
def Cauch(mean,gam):
	ran1 = n.random.uniform(0,1)
	cauch = mean + gam*n.tan(n.pi*(ran1-0.5))
	return cauch

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



def E_hf(J,I,F,A,B):#energy of each hyperfine level in eV. From allen leary's thesis
	Kk = K(F,I,J)
	E_hfine1 = (hbar*2.0*n.pi)*(0.5*(A*10**6)*Kk)#here A is converted to Hz, same for B on the next line
	E_hfine2 = (hbar*2.0*n.pi)*((B*10**6)/4.0)*(1.5*Kk*(Kk+1.0)-2.0*I*J*(I+1.0)*(J+1.0))/(I*J*(2.0*I-1.0)*(2.0*J-1.0))
	ch = m.isnan(E_hfine2)
	if m.isnan(E_hfine2):
		E_hfine2 = 0
	E = E_hfine1+E_hfine2
	return E

def Rad(N_u,J_u,N_l,J_l):#this function returns the radial matrix element for two states in units of a_0
	lw = (N_l,J_l-0.5)
	for ch in n.arange(0,14,1):
		rads = rad_mat[ch]
		if lw in [rads[0]]:
			break
	if N_u == 2:
		mat_el = rads[1]
	if N_u == 3:
		mat_el = rads[2]
	if N_u == 4:
		mat_el = rads[3]
	if N_u == 5:
		mat_el = rads[4]
	return mat_el


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
				trans.append(u-I)
				trans.append(l-I)
				amp = (2.0*l+1.0)*(2.0*u+1.0)*(wigner_6j(u,l,1,J_l,J_u,I))**2
				trans.append(amp)
				ntrans = ntrans + 1
			if ch == 0:
				if u != 0:
					trans.append(u)
					trans.append(l)
					trans.append(n.abs(u-I))
					trans.append(n.abs(l-I))
					amp = (2.0*l+1.0)*(2.0*u+1.0)*(wigner_6j(u,l,1,J_l,J_u,I))**2
					trans.append(amp)
					ntrans = ntrans + 1
	Ftemp = n.array(trans)
	F = Ftemp.reshape((ntrans,5))
	F[:,4] = F[:,4]/n.sum(F[:,4])
	return F

def mu(N_u,J_u,F_u,F_uz,N_l,J_l,F_l,F_lz,I,q):#compute the dipole moment. from eq. 4.33 pg 55 in Laser cooling and Trapping (Metcalf)
	rad_el = Rad(N_u,J_u,N_l,J_l)*a_0 #q is the photon polarization q = +-1 or 0
	exp = c_e*(-1.0)**(1+(J_u-0.5)+0.5+J_l+J_u+I-F_uz)
	sqrt = n.sqrt((2.0*J_u+1.0)*(2.0*J_l+1.0)*(2.0*F_u+1.0)*(2.0*F_l+1.0))
	wig6 = wigner_6j((J_u-0.5),J_u,0.5,J_l,(J_l-0.5),1)*wigner_6j(J_u,F_u,I,F_l,J_l,1)
	wig3 = wigner_3j(F_l,1,F_u,F_lz,q,-F_uz)
	mu = rad_el*exp*sqrt*wig6*wig3
	return mu
def A_calc(J_u,F_u,F_uz,J_l,F_l,F_lz,I,q):#compute the dipole moment. from eq. 4.33 pg 55 in Laser cooling and Trapping (Metcalf)
	#q is the photon polarization q = +-1 or 0
	exp = c_e*(-1.0)**(1+(J_u-0.5)+0.5+J_l+J_u+I-F_uz)
	sqrt = n.sqrt((2.0*J_u+1.0)*(2.0*J_l+1.0)*(2.0*F_u+1.0)*(2.0*F_l+1.0))
	wig6 = wigner_6j((J_u-0.5),J_u,0.5,J_l,(J_l-0.5),1)*wigner_6j(J_u,F_u,I,F_l,J_l,1)
	wig3 = float(wigner_3j(F_l,1,F_u,F_lz,q,-F_uz))
	A_calc = exp*sqrt*wig6*wig3
	return [A_calc,exp,sqrt,wig6,wig3]
	
def Dop(v,wn):
	beta = v/c #beta factor
	fact = ((1.0+beta)/(1.0-beta))**(0.5)
	obs = fact*wn
	return obs

def Scatter(gamma,s_0,delta):#scattering rate of photons from laser metcalf pg.25
	gamma_prime = gamma*n.sqrt(1.0+s_0)
	scattering_rate = (s_0/(1.0+s_0))*(gamma/2.0)/(1.0+(2.0*delta/gamma_prime)**2)
	return scattering_rate

def Lorentz(xdata,cent,gamma):
	Lor = (1.0/n.pi)*(gamma/((xdata-cent)**2+gamma**2))
	return Lor

def gaussian(xdata,cent,sig):
	Gauss = n.exp(-(xdata-cent)**2/(2*sig**2))/(sig*n.sqrt(2.0*n.pi))
	return Gauss


def realvoigt(xdata,cent,f_g,f_l):
	z = (xdata-cent+1j*f_l)/(f_g*n.sqrt(2*n.pi))
	wof = wofz(z)
	voigt = wof.real/(f_g*n.sqrt(2*n.pi))
	return voigt
	
	
def pseudovoigt(xdata,cent,f_g,f_l):
	f = (f_g**5+2.69269*f_g**4*f_l+2.42843*f_g**3*f_l**2+4.47163*f_g**2*f_l**3+0.07842*f_g*f_l**4+f_l**5)**(1.0/5)
	#eta =1.0 - n.abs(1-(1.36606*(f_l/f)-0.47719*(f_l/f)**2+0.11116*(f_l/f)**3))
	#eta=1
	eta = 1.36606*(f_l/f)-0.47719*(f_l/f)**2+0.11116*(f_l/f)**3
	pseudovoigt = (eta*Lorentz(xdata,cent,f_l/2.0)+(1.0-eta)*gaussian(xdata,cent,f_g/2))
	return pseudovoigt,eta

def voigt(xdata,amp,cent,sig,ep):#define the voigt peak. If ep is less than 0 or greater than 1, mess up the profile
    x = xdata
    C = cent
    S = sig
    A = amp
    E = ep
    Sg = S/n.sqrt(2*n.log(2))
    vmodel = (A*(1.0-E)/(Sg*n.sqrt(2.0*n.pi)))*n.exp(-(x-C)**2/(2.0*Sg**2))+ ((E*A)/(n.pi))*(S)/((x-C)**2+S**2)
    if 0>E or 1<E:
    	vmodel = vmodel + 50000
    return vmodel



#-----------------------Inputs-------------------------------------


v_c = raw_input('Use the same parameter vector?(y/n)')#makes it easier to rerun things

#v_c = 'y'
if v_c in ['n']:
	mss = float(raw_input('Isotope mass (AMU):'))#get mass of isotope
	Z = float(raw_input('Z (# of protons):'))#get charge of nucleus
	J_l = float(raw_input('Lower J state:'))#get lower J state of transition
	J_u = float(raw_input('Upper J state:'))#get upper J state of transition
	I = float(raw_input('Nuclear spin:'))#get nuclear spin of isotope
	gamma_e = float(raw_input('Experimental decay rate of fine structure transition (Hz):'))
	A_u = float(raw_input('Upper A coefficient (MHz):'))#hyperfine coefficients
	A_l = float(raw_input('Lower A coefficient (MHz):'))
	B_u = float(raw_input('Upper B coefficient (MHz):'))
	B_l = float(raw_input('Lower B Coefficiet (MHz):')) 
	T = float(raw_input('Beam Temperature (K):'))#get beam temperature
	l_wn = float(raw_input('Laser wavenumber (cm^(-1)):'))#get laser wavelength
	l_p = float(raw_input('Laser Power (mW):'))#get laser power
	wn_fs = float(raw_input('Fine Structure transition wavenumber (cm^(-1)):'))
	d = float(raw_input('Distance between CEC and LCR (m):'))#distance over which optical pumping is problematic
	FILENAME = str(raw_input('Experimental file (if none, then write n)'))
	vec = [mss,Z,J_l,J_u,I,gamma_e,A_u,A_l,B_u,B_l,T,l_wn,l_p,wn_fs,d,FILENAME]


if v_c in ['y']:
	mss = float(vec[0])
	Z = float(vec[1])
	J_l = float(vec[2])
	J_u = float(vec[3])
	I = float(vec[4])
	gamma_e = float(vec[5])
	A_u = float(vec[6])
	A_l = float(vec[7])
	B_u = float(vec[8])
	B_l = float(vec[9])
	T = float(vec[10])
	l_wn = float(vec[11])
	l_p = float(vec[12])
	wn_fs = float(vec[13])
	d = float(vec[14])
	FILENAME = vec[15]

m_is = amu*int(mss) #mass of the isotope
m_is_kg = amu_kg*int(mss)
En_fs = hbar*n.pi*2*c*10**2*(wn_fs)#energy in eV of fine structure
prefact = 3.0*n.pi*e_0*hbar*c**3*gamma_e/c_e**2

#Determine lifetimes and widths of each transition

F_t = F(I,J_u,J_l)#get the transitions and transition amplitudes
n_tr = F_t.shape[0]#number of transitions
n_par = F_t.shape[1]#number of parameters in F_t
F_n = n.zeros((n_tr,15))#create new matrix with extra columns of Photon energy, Energy width, lifetime, saturation intensity, saturation parameter 
F_n[0:n_tr,0:n_par] = F_t #put in original F_t

#the next section will be filling F_n with the relevant quantities, up until A


A_vec = n.zeros((n_tr))
omeg_vec = n.zeros((n_tr))
#count = 0
for k in n.arange(0,n_tr,1):
	A_calc_int=0
	A_calc_tot = 0
	#first, photon energy
	f_u = F_n[k,0]#get upper F state
	f_l = F_n[k,1]#get lower F State
	f_u_en =  E_hf(J_u,I,f_u,A_u,B_u) #get energy of upper state knowing the fine structure energy
	f_l_en =  E_hf(J_l,I,f_l,A_l,B_l) #get energy of lower state knowing the fine structure energy
	del_en = f_u_en-f_l_en + En_fs
	F_n[k,4] = del_en #get energy difference in eV and store it 
	omeg = del_en/hbar #angular frequency in Hz
	
	#next, get mean lifetime of decay. To do this, we need to compute A for each (f_uz,f_lz,q) This is spontaneous emission
	A_calc_int = 0
	A_calc_tot = 0
	for fuz in n.arange(-f_u,f_u+1,1):#for each projection on z axis of top F state
		for flz in n.arange(-f_l,f_l+1):#for each projection on z axis of lower F state
			for q in [-1.0,0,1.0]:
				A_calc_int = n.abs(A_calc(J_u,f_u,fuz,J_l,f_l,flz,I,q)[0])
				A_calc_tot = A_calc_tot+A_calc_int
	A_vec[k] = A_calc_tot
	omeg_vec[k] = omeg

R_exp = (prefact*n.sum(1.0/(A_vec**2*omeg_vec**3.0)))**(0.5)

for k in n.arange(0,n_tr,1):
	
	mutot = R_exp*c_e*A_vec[k]
	dec_r = (omeg_vec[k]**3)*(mutot**2)/((3.0*n.pi)*e_0*hbar*c**3)/10#decay rate in Hz
	tau = 1.0/dec_r
	spread = hbar/(2*tau) 
	F_n[k,5] = spread #input energy spread in eV
	F_n[k,6] = tau #input lifetime in seconds
	F_n[k,7] = F_t[k,4]
	wl_fs = 2.0*n.pi*hbar*c/En_fs
	wl = 2.0*n.pi*hbar*c/del_en #transition wavelength in meters
	I_s =  n.pi*h_js*c/(3.0*(wl**3.0)*tau) #compute saturation intensity pg. 25 in Metcalf, in W
	F_n[k,8] = I_s
	F_n[k,9] = (l_p)/(I_s/10.0) #compute on resonance saturation parameter pg. 25 in metcalf
	#compute the velocity of the atoms on resonance
	e_trans = F_n[k,4]#transition energy
	f_trans = e_trans/(hbar*2.0*n.pi)#transition frequency
	f_init = c*(l_wn*100)#initial frequency in hz
	v_trans = c*((f_trans**2-f_init**2)/(f_trans**2+f_init**2))
	F_n[k,10] = n.abs(v_trans) #speed of the atom necessary for resonance

#next, evaluate the likelihood that an atom in gs G will come back to that state after being excited
F_n_shape=F_n.shape
for k in n.arange(0,n_tr,1):
	F_n_temp = []
	trans = 0
	gs = F_n[k,1]#get the ground state
	es = F_n[k,0]#get the excited state
	for kk in n.arange(0,n_tr,1):
		if F_n[kk,0] == es:
			F_n_temp.append(F_n[kk,:])#get all the transitions with es as excited state
			trans = trans+1
	F_n_temp = n.reshape(F_n_temp,(trans,F_n_shape[1]))#reshape the array
	decay_rates = 1.0/F_n_temp[:,6]
	tot_decay_rate = n.sum(decay_rates)
	decay_rates_norm = decay_rates/tot_decay_rate
	for kkk in n.arange(0,len(decay_rates),1):#find the one transition that allows for the same gs
		if F_n_temp[kkk,1]==gs:
			prob_of_og=decay_rates_norm[kkk]#likelihood of returning to og ground state
			F_n[k,11] = decay_rates_norm[kkk]#store likelihood of returning to og ground state
			scatter_res = Scatter(1.0/F_n[k,6],F_n[k,9],0)#get scattering rate on resonance
			F_n[k,12]=scatter_res#store scattering rate on resonance
			tot_time = F_n[k,6]+1.0/scatter_res#get total transition time = scattering time + lifetime
			num_trans = m.floor((d/F_n[k,10])/tot_time)#total number of transitions rounded down
			F_n[k,13] = num_trans#store number of transitions
			prob = prob_of_og**(num_trans)#likelihood of finding og gs when you enter LCR
			F_n[k,14] = prob#store the above prob


#now start plotting what a regular spectrum should look like


xdata = n.arange(n.min(F_n[:,4])*0.999999,n.max(F_n[:,4])*1.000001,1e-8)

spec_sum = n.zeros((len(xdata),))
spec_int = 0
amp=1
#print vec[12]
#print vec[10]
for spec in n.arange(0,n_tr,1):
	gamma_stim_cal = (1.0/((1.0/(F_n[spec,6])))*n.sqrt(1.0+I/F_n[spec,8]))
	gamma_stim = 1/(1.0/((1.0/gamma_e))*n.sqrt(1.0+I/F_n[spec,8]))
	#gamma_stim_used = 1/(n.pi*2*(1/gamma_e))*hbar*2.0*n.pi*n.sqrt(1.0+I_s/F_n[spec,8])
	boltz_FWHM = F_n[spec,4] * (8.0*kB_J*T*n.log(2)/(m_is_kg*c**2))**(0.5)
	boltz_HWHM_rel = (1.0/8.0)*(2.0*kB_J*T*n.log(2)/(m_is_kg*c**2))**(0.5)
	spec_int = F_n[spec,14]*F_n[spec,7]*realvoigt(xdata,F_n[spec,4],boltz_HWHM_rel,gamma_stim_cal)#plot the peak weighed by the chances of it not being pumped
	print boltz_HWHM_rel
	print gamma_stim_cal
	#print F_n[spec,5]*n.sqrt(1.0+I/F_n[spec,8])
	#print F_n[:,14]
	spec_sum=spec_sum+spec_int


do = 1
if FILENAME in ['n']:
	do = 0
	


if do == 1:
	offset =l_wn#laser frequency in cm^-1
	lam = (1.0/offset)/100 #wl in m
	En_offset = 2.0*n.pi*hbar*c/lam#offset energy in eV
		
	data = n.loadtxt(open(FILENAME,"r"))#import the datafile with name FILENAME
	exp_x = data[:,0]
	exp_y = data[:,1]
	
	exp_x_conv = exp_x*10**6#convert to Hz
	exp_x_conv = exp_x_conv*2.0*n.pi*hbar#convert to eV
	exp_x_conv = exp_x_conv+En_offset #recenter the whole thing
	
	index = 0
	for ch in exp_y:#find the highest peak in the experimental data
		if ch == n.max(exp_y):
			max_in_exp = int(index)
			break
		index = index+1.0
	
	index = 0
	for ch in spec_sum:#find the highest peak in the experimental data
		if ch == n.max(spec_sum):
			max_in_sim = int(index)
			break
		index = index+1.0
	
	bg = 1.755167
	diff = exp_x_conv[max_in_exp] - xdata[max_in_sim]#difference in locations of highest peak
	exp_x_conv = exp_x_conv -diff#line up the two highest peaks
	offset_exp_y = exp_y-bg
	spec_sum_exp = n.zeros((len(exp_x_conv),))
	spec_int = 0
	amp=1
	for spec in n.arange(0,n_tr,1):#remake the spectrum, ready for comparison
		gamma_stim_cal = (1.0/((1.0/(F_n[spec,6])))*n.sqrt(1.0+I/F_n[spec,8]))
		gamma_stim = 1/(1.0/((1.0/gamma_e))*n.sqrt(1.0+I/F_n[spec,8]))
		#gamma_stim_used = 1/(n.pi*2*(1/gamma_e))*hbar*2.0*n.pi*n.sqrt(1.0+I_s/F_n[spec,8])
		boltz_FWHM = F_n[spec,4] * (8.0*kB_J*T*n.log(2)/(m_is_kg*c**2))**(0.5)
		boltz_HWHM_rel = (1.0/8)*(2.0*kB_J*T*n.log(2)/(m_is_kg*c**2))**(0.5)
		spec_int = F_n[spec,14]*F_n[spec,7]*realvoigt(exp_x_conv,F_n[spec,4],boltz_HWHM_rel,gamma_stim_cal)#plot the peak weighed by the chances of it not being pumped
		spec_sum_exp=spec_sum_exp+spec_int
	
	
	spec_sum_norm = spec_sum_exp/n.max(spec_sum_exp)
	mult = n.max(offset_exp_y)
	spec_sum_plot = spec_sum_norm*mult+bg
	uncert = n.sqrt(exp_y)
	index = 0#correct for no counts, which will be important
	
	for ch in uncert:
		if ch == 0:
			uncert[index] = 1
		index = index+1
	#next, remake the spectrum with exp_x_conv
	chi_2_red = 1.0/len(exp_x_conv)*n.sum((spec_sum_plot-bg-offset_exp_y)**2/uncert**2)
	plt.clf()
	plt.errorbar(exp_x_conv,exp_y,uncert,fmt='k.',label=('Ga-69 Exp.'))
	plt.plot(exp_x_conv,spec_sum_plot+bg,'-r',label=('Ga-69 Sim.'),)
	plt.legend()
	plt.xlabel('Energy (eV)')
	plt.ylabel('Photon Counts')
	plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
	plt.minorticks_on()
	plt.text(3.0e-5+2.97093,40,'$\chi^2_{red}$=%f'%(chi_2_red))
	plt.text(3.0e-5+2.97093,50,'$Power = %f$\mu$W'%(vec[12]/10e-3))
	plt.ylim((-1,n.max(spec_sum_plot)))
	plt.savefig('Ga-69-vs-sim.png')
	print chi_2_red	
	
	




