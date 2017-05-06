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
					trans.append(u-I)
					trans.append(l-I)
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
	vec = [mss,Z,J_l,J_u,N_l,N_u,I,A_u,A_l,B_u,B_l,T,l_wn,l_p,wn_fs,t_v,dac_step,d,d_LCR]
	#print '[%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f]'%(mss,Z,J_l,J_u,N_l,N_u,I,A_u,A_l,B_u,A_l,T,l_f,v_i,v_inc,t_v,d)


if v_c in ['y']:
	mss = float(vec[0])
	Z = float(vec[1])
	J_l = float(vec[2])
	J_u = float(vec[3])
	N_l = float(vec[4])
	N_u = float(vec[5])
	I = float(vec[6])
	A_u = float(vec[7])
	A_l = float(vec[8])
	B_u = float(vec[9])
	B_l = float(vec[10])
	T = float(vec[11])
	l_wn = float(vec[12])
	l_p = float(vec[13])
	wn_fs = float(vec[14])
	t_v = float(vec[15])
	dac_step = float(vec[16])
	d = float(vec[17])
	d_LCR = float(vec[18])


#fn_st = raw_input('Simulate the fine structure energy?(y/n)')
#if fn_st in ['n']:
#	wn_fs = float(raw_input('Wavenumber of fine structure transition (cm^(-1)):'))
#	En_fs = hbar*n.pi*2*c*10**2*(wn_fs)#energy in eV of fine structure


#--------------Actual Simulation begins here--------------------

#preliminaries
m_is = amu*int(mss) #mass of the isotope
En_fs = hbar*n.pi*2*c*10**2*(wn_fs)#energy in MeV of fine structure

#Determine lifetimes and widths of each transition

F_t = F(I,J_u,J_l)#get the transitions and transition amplitudes
n_tr = F_t.shape[0]#number of transitions
n_par = F_t.shape[1]#number of parameters in F_t
F_n = n.zeros((n_tr,n_par+5))#create new matrix with extra columns of Photon energy, Energy width, lifetime, saturation intensity, saturation parameter 
F_n[0:n_tr,0:n_par] = F_t #put in original F_t

#the next section will be filling F_n with the relevant quantities

for k in n.arange(0,n_tr,1):
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
	
	#next, get mean lifetime of decay. To do this, we need to compute mu for each (f_uz,f_lz,q) This is spontaneous emission
	mutot_sqrd = 0
	mutot = 0
	for fuz in n.arange(-f_u,f_u+1,1):#for each projection on z axis of top F state
		for flz in n.arange(-f_l,f_l+1):#for each projection on z axis of lower F state
			for q in [-1.0,0,1.0]:
				mu_int = mu(N_u,J_u,f_u,fuz,N_l,J_l,f_l,flz,I,q)
				mutot = mutot+mu_int
				mutot_sqrd = mu_int**2 + mutot_sqrd
				
	
	dec_r = (omeg**3)*(mutot_sqrd)/((3.0*n.pi)*e_0*hbar*c**3)#decay rate in Hz
	tau = 1.0/(0.945*10**8) #lifetime in seconds...average. Got this from
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
#vs = [] #store the velocities of the atoms
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

dacs = n.arange((beam_en_min*0.9995),(beam_en_max*1.0005+dac_step),(dac_step))#beam energy steps in eV
l = float(len(dacs))
ph_detected = 0.0
prog_dac = 0

#-----not right--------
#r_mean = N_l**2*a_0/Z#mean radius in m
#crss_area = n.pi*r_mean**2 #cross sectional area of atom
#l_p_W = float(l_p)/10**3 #laser power in Watts
#l_en = h_js*c/(1.0/(100.0*l_wn))#laser energy in joules
#photon_rate = l_p_W/l_en

#this is the proper scattering rate



print 'Starting run...'
ti.sleep(1)
for jk in dacs:#step atoms through dac steps
	ph_dac=0 #photons detected at each dac step
	#if v_i+jk <= 0:
	#	v_beam = -c*n.sqrt(2.*(n.abs(v_i+jk)/(m_is)))#get the speed of the beam in m/s after the acceleration due to the doppler voltage
	#if v_i+jk > 0:
	#	v_beam = c*n.sqrt(2.*(n.abs(v_i+jk)/(m_is)))
	v_beam = c*n.sqrt(2.*(n.abs(jk)/(m_is)))

	for i in n.arange(0,t_v,1):#this is the loop that will run each atom through the optical pumping region and the LCR
		#v_therm = Boltz(m_is,T)#give the atom a shift in on axis velocity according to boltzmann distribution
		v_therm = 0
		vs_therm.append(v_therm)#store the shift
		v_a = v_beam + v_therm #final velocity of the atom
		vs.append(v_a)
		
		
		l_wnd = Dop(v_a,l_wn)#get the shifted wavenumber of the laser in cm^-1
		l_lam = 1/(100*l_wnd)#wavelength of shifted laser in m
		
		#------Ignore------#
		#int_cross = 4.0*(alpha**4)*n.sqrt(2.0)*(Z**5)*(8.0/3.0)*n.pi*(r_e**2)*((m_e)/(hbar*2.0*n.pi*c/l_lam))**(7.0/2.0)#interaction cross section for atom: http://www.physics.queensu.ca/~phys352/lect17.pdf
		#ph_sec = int_cross*photon_rate#number of photons per second
		#t_int = 1.0/(ph_sec) #time per exciting photon
		#stop ignoring-----#
		
		
		wn.append(l_wnd)
		en_l = hbar*2.0*n.pi*c*l_wnd*10**2#get the energy of the shifted laser in eV
		ens.append(en_l)#store the laser energies
		
		
		#now, we decide which F state the atom enters with
		F_g = GSP(I,J_l)#atom is now in a particular ground state
		F_ex = 0
		gs.append(F_g) #store all the ground states of the incoming atoms just to make sure everything is ok
		 #skip atom or not. Happens if there are no available initial transitions. no skip = 0, skip = 1
		pos = 0 #position of the ion
		trans_en = F_n[:,4:6] #transition energy and widths
		
		while pos < d:#Here the atom will be excited and deexcited in the distance between the CEC and the LCR
			skip = 1#going to skip unless we get a transition
			
			for trans in n.arange(0,n_tr,1):#first see if any transitions are available
				
				if trans_en[trans,0]-trans_en[trans,1]<=en_l<=trans_en[trans,0]+trans_en[trans,1] and F_n[trans,1] == F_g: #if the transition is available, transition
					delta = n.abs(en_l-trans_en[trans,0])/(hbar*2.0*n.pi)
					scatter_rate = Scatter(1.0/F_n[trans,6],F_n[trans,9],delta)
					t_int = 1.0/scatter_rate
					pos = pos+t_int*v_a#advance the atom by the time it takes to interact with a photon
					F_ex = F_n[trans,0] #set the state to the upper state of the transition
					#print 'Excite'
					#Next we de-excite the atom
					de_ex_mat = []
					n_trans_de_ex = 0
					
					for de_ex in n.arange(0,n_tr,1):#cycle through the transitions, looking for those with the allowed lower states
						if F_n[de_ex,0] == F_ex: #select the transitions with F_ex as the upper state
							n_trans_de_ex = n_trans_de_ex + 1 #increase number of available transitions
							de_ex_mat.append(F_n[de_ex,:])#store that potential transition
							
							
					de_ex_mat = n.reshape(de_ex_mat,(n_trans_de_ex,10))#reshape the deexcitation matrix
					tot_lik = n.sum(de_ex_mat[:,7]) #use lifetimes to decide total likelihood of decay
					
					likelihoods = (de_ex_mat[:,7])/tot_lik #likelihood of this decay happening
					r_num = r.uniform(0,1)#generate a random number between 0 and 1
					cumulative = 0 #cumulative probability
					
					#choose the transition that will happen
					for de_ex_choice in n.arange(0,n_trans_de_ex,1):
						prob = likelihoods[de_ex_choice]
						if cumulative <= r_num <= prob+cumulative:
							transition = de_ex_mat[de_ex_choice,:]#if the random number falls in the right range, pick that transtion
							#print 'Deexcite'
							break
						cumulative = cumulative + prob#if the random number doesnt fall in the right range, change range and check again
					Lifetime = transition[6] #get the lifetime in seconds
					pos = pos + Lifetime*v_a #advance the atom by the lifetime times the speed
					F_g = transition[1] #set the new ground state
					skip = 0
					
			if skip == 1:
				pos = d
			#----------------------------------------IN THE LCR--------------
		gs_LCR.append(F_g)#store the ground state the atom comes in as
		while d <= pos < d+d_LCR:#Here the atom will be excited and deexcited in the LCR
			skip = 1
			for trans in n.arange(0,n_tr,1):#first see if any transitions are available
			
				if trans_en[trans,0]-trans_en[trans,1]<=en_l<=trans_en[trans,0]+trans_en[trans,1] and F_n[trans,1] == F_g: #if the transition is available, transition
					#print 'TRANSITION IN LCR'
					delta = n.abs(en_l-trans_en[trans,0])/(hbar*2.0*n.pi)
					scatter_rate = Scatter(1.0/F_n[trans,6],F_n[trans,9],delta)
					t_int = 1.0/scatter_rate
					pos = pos+t_int*v_a#advance the atom by the time it takes to interact with a photon
					F_ex = F_n[trans,0] #set the state to the upper state of the transition
					#Next we de-excite the atom
					de_ex_mat = []
					n_trans_de_ex = 0
					for de_ex in n.arange(0,n_tr,1):#cycle through the transitions, looking for those with the allowed lower states
						if F_n[de_ex,0] == F_ex: #select the transitions with F_ex as the upper state
							n_trans_de_ex = n_trans_de_ex + 1 #increase number of available transitions
							de_ex_mat.append(F_n[de_ex,:])#store that potential transition
					
					de_ex_mat = n.reshape(de_ex_mat,(n_trans_de_ex,10))#reshape the deexcitation matrix
					tot_lik = n.sum(de_ex_mat[:,7]) #use lifetimes to decide total likelihood of decay
					
					likelihoods = (de_ex_mat[:,7])/tot_lik #likelihood of this decay happening
					r_num = r.uniform(0,1)#generate a random number between 0 and 1
					cumulative = 0 #cumulative probability
				
					#choose the transition that will happen
					for de_ex_choice in n.arange(0,n_trans_de_ex,1):
						prob = likelihoods[de_ex_choice]
						if cumulative <= r_num <= prob+cumulative:
							transition = de_ex_mat[de_ex_choice,:]#if the random number falls in the right range, pick that transtion
							break
						cumulative = cumulative + prob#if the random number doesnt fall in the right range, change range and check again
					
					Lifetime = transition[6] #get the lifetime in seconds
					pos = pos + Lifetime*v_a #advance the atom by the lifetime times the speed
					photons.append(transition[4]+Cauch(0,transition[5]))
					ph_dac = ph_dac+1
					F_g = transition[1] #set the new ground state
					skip = 0
			if skip == 1:
				pos = d+d_LCR

	prog = n.round(float(prog_dac/l),6)
	print '{} {}\r'.format(prog,ph_dac)#show progress
	prog_dac = prog_dac+1.0		
	cts_dac.append(ph_dac)	

with open('Output_new.txt','w') as txt_file:
	for i in cts_dac:
		txt_file.write('%f\n' % (i))

with open('Params_new.txt','w') as txt_file:
	for param in vec:
		txt_file.write('%f\n' % (param))

zer0 = 0
one = 0
two = 0
three = 0
for hist in gs_LCR:
	if hist == 0:
		zer0 = zer0+1
	
	if hist == 1:
		one = one+1
	
	if hist == 2:
		two = two+1
	
	if hist == 3:
		three = three+1

gs_pumped = n.array([zer0,one,two,three])

zer0 = 0
one = 0
two = 0
three = 0
for hist in gs:
	if hist == 0:
		zer0 = zer0+1
	
	if hist == 1:
		one = one+1
	
	if hist == 2:
		two = two+1
	
	if hist == 3:
		three = three+1

gs = n.array([zer0,one,two,three])

with open('gs_pumped.txt','w') as txt_file:
	for gs_pump in gs_pumped:
		txt_file.write('%f\n' % (gs_pump))

with open('photons.txt','w') as txt_file:
	for phot in photons:
		txt_file.write('%s\n' % (str(Decimal(phot))))

ph_hist = n.histogram(photons, bins = 1000)