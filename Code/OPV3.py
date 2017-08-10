import numpy as n
import math as m
import pyspec as p
import time as ti
from decimal import *
import random as r
import numpy.linalg as lin
from sympy.physics.wigner import wigner_6j, wigner_3j
import matplotlib.pyplot as plt


#----------------------INTRO---------------------------


#------------------Preamble------------------------------------
# Important constants

kB = 8.6173324 * 10**(-5) #eV/K boltzmanns constant
amu = 931.4940954*10**6 #eV/c^2 rest energy of one atomic mass unit
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

#define the radial transition matrix elements in units of a_o
#obtained from Laser cooling and trapping (Metcalf) pg.52
rad_mat = [[(1,0),1.2902,0.5166,0.3045,0.2087],
	[(2,0),-5.1961,3.0648,1.2822,0.7739],
	[(2,1),'NA',4.7479,1.7097,0.9750],
	[(3,0),0.9384,-12.7279,5.4693,2.25957],
	[(3,1),'NA',-10.0623,7.5654,2.9683],
	[(3,2),'NA','NA',10.2303,3.3186],
	[(4,0),0.3823,2.4435,-23.2379,8.5178],
	[(4,1),'NA',1.3022,-20.7846,11.0389],
	[(4,2),'NA','NA',-15.8745,14.0652],
	[(4,3),'NA','NA','NA',17.7206],
	[(5,0),0.2280,0.9696,4.6002,-36.7423],
	[(5,1),'NA',0.4827,3.0453,-34.3693],
	[(5,2),'NA','NA',1.6613,-30.0000],
	[(5,3),'NA','NA','NA',-22.5000]]






#-------------Definitions--------------------#


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


#def E_f(N,J,Z):#Energy of the fine structure levels in eV. From Pg. 408 in Townsend QM 
#	sqrd = Z*alpha/(N-(J+0.5)+((J+0.5)**2-(Z*alpha)**2)**(0.5))
#	E = m_e*((1.0+sqrd**2)**(-0.5)-1)
#	return E




def E_hf(J,I,F,A,B):#energy of each hyperfine level in MeV. From allen leary's thesis
	Kk = K(F,I,J)
	E_hfine1 = (hbar/(2.0*n.pi))*(0.5*(A*10**6)*Kk)#here A is converted to Hz, same for B on the next line
	E_hfine2 = (hbar/(2.0*n.pi))*((B*10**6)/4.0)*(1.5*Kk*(Kk+1.0)-2.0*I*J*(I+1.0)*(J+1.0))/(I*J*(2.0*I-1.0)*(2.0*J-1.0))
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
#-----------------------Inputs-------------------------------------


v_c = raw_input('Use the same parameter vector?(y/n)')#makes it easier to rerun things

if v_c in ['n']:
	mss = float(raw_input('Isotope mass (AMU):'))#get mass of isotope
	Z = float(raw_input('Z (# of protons):'))#get charge of nucleus
	J_l = float(raw_input('Lower J state:'))#get lower J state of transition
	J_u = float(raw_input('Upper J state:'))#get upper J state of transition
	N_l = float(raw_input('Lower N: '))#get lower energy state
	N_u = float(raw_input('Upper N: '))#get upper energy state
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
	t_v = float(raw_input('Number of atoms to simulate per DAC increment:'))#number of atoms that will be simulated
	dac_step = float(raw_input('Beam Energy Step (eV):'))#voltage step
	d = float(raw_input('Distance between CEC and LCR (m):'))#distance over which optical pumping is problematic
	d_LCR = float(raw_input('Length of LCR (m):'))#length of the LCR
	vec = [mss,Z,J_l,J_u,N_l,N_u,I,gamma_e,A_u,A_l,B_u,B_l,T,l_wn,l_p,wn_fs,t_v,dac_step,d,d_LCR]
	#print '[%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f]'%(mss,Z,J_l,J_u,N_l,N_u,I,A_u,A_l,B_u,A_l,T,l_f,v_i,v_inc,t_v,d)


if v_c in ['y']:
	mss = float(vec[0])
	Z = float(vec[1])
	J_l = float(vec[2])
	J_u = float(vec[3])
	N_l = float(vec[4])
	N_u = float(vec[5])
	I = float(vec[6])
	gamma_e = float(vec[7])
	A_u = float(vec[8])
	A_l = float(vec[9])
	B_u = float(vec[10])
	B_l = float(vec[11])
	T = float(vec[12])
	l_wn = float(vec[13])
	l_p = float(vec[14])
	wn_fs = float(vec[15])
	t_v = float(vec[16])
	dac_step = float(vec[17])
	d = float(vec[18])
	d_LCR = float(vec[19])


#fn_st = raw_input('Simulate the fine structure energy?(y/n)')
#if fn_st in ['n']:
#	wn_fs = float(raw_input('Wavenumber of fine structure transition (cm^(-1)):'))
#	En_fs = hbar*n.pi*2*c*10**2*(wn_fs)#energy in eV of fine structure


#--------------Actual Simulation begins here--------------------

#preliminaries
m_is = amu*int(mss) #mass of the isotope
En_fs = hbar*n.pi*2*c*10**2*(wn_fs)#energy in MeV of fine structure
prefact = 3.0*n.pi*e_0*hbar*c**3*gamma_e/c_e**2

#Determine lifetimes and widths of each transition

F_t = F(I,J_u,J_l)#get the transitions and transition amplitudes
n_tr = F_t.shape[0]#number of transitions
n_par = F_t.shape[1]#number of parameters in F_t
F_n = n.zeros((n_tr,n_par+5))#create new matrix with extra columns of Photon energy, Energy width, lifetime, saturation intensity, saturation parameter 
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
	#print f_u_en
	f_l_en =  E_hf(J_l,I,f_l,A_l,B_l) #get energy of lower state knowing the fine structure energy
	#print f_l_en
	del_en = f_u_en-f_l_en + En_fs
	F_n[k,4] = del_en #get energy difference in eV and store it 
	omeg = del_en/hbar #angular frequency in Hz
	
	#next, get mean lifetime of decay. To do this, we need to compute A for each (f_uz,f_lz,q) This is spontaneous emission
	A_calc_int = 0
	A_calc_tot = 0
	for fuz in n.arange(-f_u,f_u+1,1):#for each projection on z axis of top F state
		for flz in n.arange(-f_l,f_l+1):#for each projection on z axis of lower F state
			for q in [-1.0,0,1.0]:
				#mu_int = mu(N_u,J_u,f_u,fuz,N_l,J_l,f_l,flz,I,q)
				A_calc_int = n.abs(A_calc(J_u,f_u,fuz,J_l,f_l,flz,I,q)[0])
				A_calc_tot = A_calc_tot+A_calc_int
				#if A_calc_int == 0:
				#	print J_u, f_u, fuz, J_l, f_l, flz, I, q
				#	print A_calc(J_u,f_u,fuz,J_l,f_l,flz,I,q)
				#	count = count+1
				#	print count
				#mutot = mutot+mu_int
				#mutot_sqrd = mu_int**2 + mutot_sqrd
	A_vec[k] = A_calc_tot
	omeg_vec[k] = omeg

R_exp = (prefact*n.sum(1.0/(A_vec**2*omeg_vec**3.0)))**(0.5)

for k in n.arange(0,n_tr,1):
	
	mutot = R_exp*c_e*A_vec[k]
	dec_r = (omeg_vec[k]**3)*(mutot**2)/((3.0*n.pi)*e_0*hbar*c**3)/10#decay rate in Hz
	tau = 1.0/dec_r
	#tau = 1.0/(0.945*10**8) #lifetime in seconds...average. Got this from
	spread = hbar/(2*tau) 
	F_n[k,5] = spread #input energy spread in eV
	#F_n[k,5] = 1e-9
	F_n[k,6] = tau #input lifetime in seconds
	F_n[k,7] = F_t[k,4]
	wl = 2.0*n.pi*hbar*c/del_en #transition wavelength in meters
	I_s =  n.pi*h_js*c/(3*(wl**3)*tau) #compute saturation intensity pg. 25 in Metcalf, in W
	F_n[k,8] = I_s
	F_n[k,9] = (l_p/(1000.0))/I_s #compute on resonance saturation parameter pg. 25 in metcalf
	

e_min = min(F_n[:,4])#lowest energy transition in ev
e_max = max(F_n[:,4])#highest energy transition in ev
f_min = e_min/(hbar*2.0*n.pi)#lowest frequency in hz
f_max = e_max/(hbar*2.0*n.pi)#highest frequency in hz
f_init = c*(100*l_wn)#initial frequency in hz


v_min = c*((f_min**2-f_init**2)/(f_min**2+f_init**2))#minimum velocity of atoms in m/s
v_max = c*((f_max**2-f_init**2)/(f_max**2+f_init**2))#maximum velocity of atoms in m/s

beam_en_min = 0.5*m_is*(v_min/c)**2#min beam energy in eV
beam_en_max = 0.5*m_is*(v_max/c)**2#max beam energy in eV
#ok now we know everything about each transition

#Next we start figuring out what the beam looks like:

#v_beam_init = c*n.sqrt(2.*(v_i/m_is))#get the speed of the beam in m/s
vs = [] #store the velocities of the atoms
vs_therm=[]
gs = [] #store all the ground states
gs_LCR=[]#store the states that the atoms come into the LCR in
photons = [] #store the photons measured by the PMT
cts_dac = [] #count the photons at each dac increment 
ens = [] #store energies
cts = [] #store counts
wn = [] #store shifted wavenumber of laser
uppers = [] #store energy limits
lowers = []
scatter_rates_store = []
scatter_times = []
deltas=[]
states_hist = []

diff = beam_en_max-beam_en_min
dacs = n.arange((beam_en_min-0.01*diff),(beam_en_max+0.01*diff+dac_step),(dac_step))#beam energy steps in eV
l = float(len(dacs))
ph_detected = 0.0
prog_dac = 0

gs_array = n.zeros((len(dacs),4,2))
ct_gs = 0
print 'Starting run...'

ti.sleep(1)

for jk in dacs:
	ph_dac = 0 #Number of photons detected at each dac increment set to 0
	v_beam = c*n.sqrt(2.*(n.abs(jk)/(m_is)))#find the velocity of the beam at this DAC increment
	
	for i in n.arange(0,t_v,1):#this loops through every atom at each dac increment
		states_inc = []
		v_therm = Boltz(m_is,T)#give the atom a shift in on axis velocity according to boltzmann distribution depending on the mass and temperature of the beam
		vs_therm.append(v_therm)#store the shift
		v_a = v_beam + v_therm #final velocity of the atom
		vs.append(v_a)
		
		l_wnd = Dop(v_a,l_wn)#get the shifted wavenumber of the laser in cm^-1
		l_lam = 1/(100*l_wnd)#wavelength of shifted laser in m
		
		wn.append(l_wnd)#store the shifted wavelength of the laser
		en_l = hbar*2.0*n.pi*c*l_wnd*10**2#get the energy of the shifted laser in eV
		ens.append(en_l)#store the laser energies
		F_g = GSP(I,J_l)#give the atom a ground F state
		states_inc.append(F_g)
		F_ex = 0 #define and initialize the excited state
		gs.append(F_g) #store all the ground states of the incoming atoms just to make sure everything is ok
		pos = 0 #set the position of the atom to the beginning of the CEC
		trans_en = F_n[:,4:5]#get the energy of the transitions as well as the spread for deexcitation
		
		while pos < d+d_LCR: #Here the atom is in the CEC-LCR space.
		#THIS IS THE EXCITATION CODE
		#First, find all available transitions given F_g
			F_trans = []
			trans_count = 0
			for f_trans in n.arange(0,n_tr,1):
				F_g_int = F_n[f_trans,1]#find the ground state of a transition
				if F_g_int == F_g:
					F_trans.append(F_n[f_trans,:])
					trans_count = trans_count+1
			F_trans = n.reshape(F_trans,(trans_count,n_par+5))#now have info on all possible transitions
			#next look at likelihood of transitions given energy of beam
			delta = n.abs(en_l-F_trans[:,4])/(2*n.pi*hbar)#difference between frequency of laser and transition energies
			for nr in delta:
				deltas.append(nr)
			#next, calculate scattering rates for each transition
			scatter_rates = n.zeros((trans_count,1))
			for scatter in n.arange(0,trans_count,1):
				scatter_calc = Scatter(1.0/F_trans[scatter,6],F_trans[scatter,9],delta[scatter])
				scatter_rates[scatter] = scatter_calc
				scatter_rates_store.append(scatter_calc)
			#next find the average scatter time
			scatter_time = 1.0/(n.sum(1.0/scatter_rates))/1000000
			scatter_times.append(scatter_time)
			pos = pos + scatter_time*v_a
			if pos > d+d_LCR:
				#print 'LEFT IR'
				break
			#print 'EXCITED' #if we've left the CEC-LCR region, move to next part of the code
			#now we choose which state to select
			scatter_rates = scatter_rates/(n.sum(scatter_rates))#normalize scatter rates into likelihoods
			rand = r.uniform(0,1)
			prob = 0
			for exc in n.arange(0,trans_count,1):
				cumulative = prob + scatter_rates[exc]
				if prob <= rand < cumulative:
					F_exc = F_trans[exc,0]
					states_inc.append(F_exc)
					break
				prob = cumulative
			#THIS IS THE DE-EXCITATION CODE
			F_trans = []
			trans_count = 0
			for f_trans in n.arange(0,n_tr,1):
				F_exc_int = F_n[f_trans,0]#find the excited state of a transition
				if F_exc_int == F_exc:
					F_trans.append(F_n[f_trans,:])
					trans_count = trans_count+1
			F_trans = n.reshape(F_trans,(trans_count,n_par+5))#now have info on all possible transitions
			lifetimes = F_trans[:,6]
			de_excite_time = 1.0/n.sum(1.0/decay_rates)#get the cumulative lifetime of the decay
			pos = pos+ de_excite_time*v_a#advance by the lifetime of the state
			#if pos > d+d_LCR:
			#	break
			decay_rates = decay_rates/(n.sum(decay_rates))#normalize the decay rates
			rand = r.uniform(0,1)
			prob = 0 
			for exc in n.arange(0,trans_count,1):
				cumulative = prob + scatter_rates[exc]
				if prob <= rand < cumulative:
					F_g = F_trans[exc,1]
					states_inc.append(F_g)
					photon_emit = F_trans[exc,4]+Cauch(0,F_trans[exc,5])#get the energy of the photon released by the decay
					if d<=pos<=d_LCR:#if the atom is in the LCR, record the photon
						photons.append(photon_emit)
						print'Photon'
					break
				prob = cumulative
		states_hist.append(states_inc)
	prog = n.round(float(prog_dac/l),6)
	print '{} {}\r'.format(prog,ph_dac)#show progress
	prog_dac = prog_dac+1.0		
	cts_dac.append(ph_dac)
				
		
		
		
		
		
		