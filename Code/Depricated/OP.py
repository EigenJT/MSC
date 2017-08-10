import numpy as n
import math as m
import pyspec as p
import time as ti
import random as r
import numpy.linalg as lin
from sympy.physics.wigner import wigner_6j, wigner_3j
import matplotlib.pyplot as plt


#----------------------INTRO---------------------------


#------------------Preamble------------------------------------
# Important constants

kB = 8.6173324 * 10**(-11) #MeV/K boltzmanns constant
amu = 931.4940954 #MeV/c^2 rest energy of one atomic mass unit
c = 299792458 #m/s speed of light
alpha = 0.0072973525664 #fine structure constant
m_e = 0.5109989461#MeV/c^2 mass of electron
c_e = 1.6021766208*10**(-19) #charge of electron in Coulombs
kB_cc = kB/(c**2)#kB/c^2 boltzmann constant per c^2
hbar = 6.582119514*10**(-22) #reduced planck constant in MeV*s
e_0 = 8.854187817*10**(-6)/(6.242*10**(18)) #vacuum perimittivity in Coulombs^2/MeV*m
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
	wn_fs = float(raw_input('Fine Structure transition wavenumber (cm^(-1)):'))
	#l_p = float(raw_input('Laser Drift (nm):'))#get laser drift
	v_i = float(raw_input('Inital Beam Energy (MeV):'))#initial beam energy
	#v_inc = float(raw_input('Acceleration Voltage increment (V)'))#get increment of acceleration voltage
	t_v = float(raw_input('Number of atoms to simulate per DAC increment:'))#number of atoms that will be simulated
	dac_init = float(raw_input('Start Voltage (MeV):'))#initial voltage
	dac_end = float(raw_input('End Voltage (MeV):'))#final voltage
	dac_step = float(raw_input('Voltage Step (MeV):'))#voltage step
	d = float(raw_input('Distance between CEC and LCR (m):'))#distance over which optical pumping is problematic
	d_LCR = float(raw_input('Length of LCR (m):'))#length of the LCR
	vec =[mss,Z,J_l,J_u,N_u,N_l,I,A_u,A_l,B_u,B_l,T,l_wn,wn_fs,v_i,t_v,dac_init,dac_end,dac_step,d,d_LCR]
	#print '[%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f]'%(mss,Z,J_l,J_u,N_l,N_u,I,A_u,A_l,B_u,A_l,T,l_f,v_i,v_inc,t_v,d)


if v_c in ['y']:
	mss = vec[0]
	Z = vec[1]
	J_l = vec[2]
	J_u = vec[3]
	N_l = vec[4]
	N_u = vec[5]
	I = vec[6]
	A_u = vec[7]
	A_l = vec[8]
	B_u = vec[9]
	B_l = vec[10]
	T = vec[11]
	l_wn = vec[12]
	wn_fs = vec[13]
	v_i = vec[14]
	#v_inc = vec[12]
	t_v = vec[15]
	dac_init = vec[16]
	dac_end = vec[17]
	dac_step = vec[18]
	d = vec[19]
	d_LCR = vec[20]

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
F_n = n.zeros((n_tr,n_par+3))#create new matrix with extra columns of Photon energy, Energy width and lifetime. 
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
	F_n[k,4] = del_en #get energy difference in MeV and store it 
	omeg = del_en/hbar #angular frequency in Hz
	
	#next, get mean lifetime of decay. To do this, we need to compute mu for each (f_uz,f_lz,q) This is spontaneous emission
	mutot_sqrd = 0
	for fuz in n.arange(-f_u,f_u+1,1):#for each projection on z axis of top F state
		for flz in n.arange(-f_l,f_l+1):#for each projection on z axis of lower F state
			for q in [-1.0,0,1.0]:
				mu_int = mu(N_u,J_u,f_u,fuz,N_l,J_l,f_l,flz,I,q)
				mutot_sqrd = mu_int**2 + mutot_sqrd
				#print mu_int/(a_0*c_e)
	
	dec_r = (omeg**3)*(mutot_sqrd)/((3.0*n.pi)*e_0*hbar*c**3)#decay rate in Hz
	tau = 1.0/dec_r #lifetime in seconds
	spread = hbar/(2*tau) 
	F_n[k,5] = spread #input energy spread in MeV
	F_n[k,6] = tau #input lifetime in seconds
	F_n[k,7] = F_t[k,4]
	


#ok now we know everything about each transition

#Next we start figuring out what the beam looks like:

v_beam_init = c*n.sqrt(2.*(v_i/m_is))#get the speed of the beam in m/s
vs = [] #store the velocities of the atoms
vs_therm=[]
gs = [] #store all the ground states
gs_LCR=[]#store the states that the atoms come into the LCR in
photons = [] #store the photons measured by the PMT
cts_dac = [] #count the photons at each dac increment 
ens = [] #store energies
cts = [] #store counts
wn = [] #store wn
uppers = [] #store energy limits
lowers = []

dacs = n.arange(dac_init,dac_end+dac_step,dac_step)
jjj = 0.0
ph_d = 0.0
for jk in dacs:#step atoms through dac steps
	l = float(len(dacs))
	prog = n.round(float(jjj/l),6)
	print '{} {}\r'.format(prog,ph_d)#show progress
	ph_d = 0.0#photons per dac step
	jjj = jjj+1.0
	
	if v_i+jk <= 0:
		v_beam = -c*n.sqrt(2.*(n.abs(v_i+jk)/(m_is)))#get the speed of the beam in m/s after the acceleration due to the doppler voltage
	if v_i+jk > 0:
		v_beam = c*n.sqrt(2.*(n.abs(v_i+jk)/(m_is)))
	
	for i in n.arange(0,t_v,1):#this is the loop that will run each atom through the optical pumping region and the LCR
		
		v_therm = Boltz(m_is,T)#give the atom a shift in on axis velocity according to boltzmann distribution
		vs_therm.append(v_therm)#store the shift
		v_a = v_beam + v_therm #final velocity of the atom
		vs.append(v_a)
		tm = t_v/v_a #time in seconds over which optical pumping is an issue
		
		
		
		l_wnd = Dop(v_a,l_wn)#get the shifted wavenumber of the laser in cm^-1
		wn.append(l_wnd)
		en_l = hbar*2.0*n.pi*c*l_wnd*10**2#get the energy of the shifted laser in MeV
		ens.append(en_l)#store the laser energies
		
		
		#now, we decide which F state the atom enters with
		F_g = GSP(I,J_l)#atom is now in a particular ground state
		gs.append(F_g) #store all the ground states of the incoming atoms just to make sure everything is ok
		al = []
		excc = 0 #excited or not. 0-> not excited, 1-> excited
		skip = 0 #skip atom or not. Happens if there are no available initial transitions. no skip = 0, skip = 1
		pos = 0 #position of the ion
		
		
		while pos <= d:#this is the loop that will excite and de-excite the atom until it reaches the LCR
			al = []
			n_trr = 0
			#---------------------------Excite--------------------------------------
			for j in n.arange(0,n_tr,1):#What excited states are available to it? 
				ch = F_n[j,1]#get the lower F state of the transition
				if ch == F_g:#find all the transitions with F_g as the lower state
					al.append(F_n[j,:])
					n_trr = n_trr + 1
			al = n.reshape(al,(n_trr,8))
			
			for j in n.arange(0,al.shape[0],1):#cycle through the available transitions
				#print (al[j,4]-al[j,5])
				#print (al[j,4]+al[j,5])

				if ((al[j,4]-al[j,5]))<= en_l <= (al[j,4]+al[j,5]):#check to see if the transition is available to the atom
					F_g = al[j,0] #if the transition is available, excite the atom to that F state
					excc = 1#excited state of atom
					#print 'Excited OPR'
				uppers.append(al[j,4]+al[j,5])
				lowers.append(al[j,4]-al[j,5])
			if excc == 0:
				skip = 1#if the atom doesn't get excited, skip it
				break#stop this while loop
			pos = pos+v_a*al[j,6]#move the atom the distance of the lifetime*v_a
			
			
			#-------------------------De-Excite-------------------------------------
			al = []
			n_trr = 0		
			for j in n.arange(0,n_tr,1):#What de-excited states are available to it?
				ch = F_n[j,0]#get the lower F state of the transition
				if ch == F_g:#find all the transitions with F_g as the upper
					al.append(F_n[j,:])#append properties of transition
					n_trr = n_trr + 1
			al = n.reshape(al,(n_trr,8))
			al[:,7] = al[:,7]/n.sum(al[:,7])#normalize transition strengths
			r_num = r.uniform(0,1) #generate a number between 0 and 1
			tot = al[0,7]
			low = 0
			for jj in n.arange(0,n_trr,1):
				if low<=r_num<=tot:
					F_g = al[jj,1]
					photon = al[jj,4] + r.uniform(-al[jj,5],al[jj,5])
				else:
					low = tot
					if jj == n_trr-1:
						tot = 1.0
					else:
						tot = tot+al[jj+1,7]
			excc = 0 #deexcite atom state
			#-------------atom is now deexcited-------------------------------			
			
			
		#--------------NOW THE ATOM IS IN THE LCR------------------#
		#this while-loop will go on until the atom gets into the LCR
		
		gs_LCR.append(F_g) #store all the ground states of the incoming atoms just to make sure everything is ok
		al = []
		
		while pos <= d+d_LCR:#this is the loop that will excite and de-excite the atom until it leaves the LCR
			if skip == 1: #quick check to see if this atom should be skipped
				break
			
			n_trr = 0
			al = []
			#---------------------------Excite--------------------------------------
			for j in n.arange(0,n_tr,1):#What excited states are available to it? 
				ch = F_n[j,1]#get the lower F state of the transition
				if ch == F_g:#find all the transitions with F_g as the lower state
					al.append(F_n[j,:])
					n_trr = n_trr + 1
			al = n.reshape(al,(n_trr,8))
			
			for j in n.arange(0,al.shape[0],1):#cycle through the available transitions
				if (al[j,4]-al[j,5])<= en_l <= (al[j,4]+al[j,5]):#check to see if the transition is available to the atom
					F_g = al[j,0] #if the transition is available, excite the atom to that F state
					excc = 1#excited state of atom
					#print 'Excited LCR'
				pos = pos+v_a*al[j,6]#move the atom the distance of the lifetime*v_a
			
			
			#-------------------------De-Excite-------------------------------------
			al = []
			n_trr = 0		
			for j in n.arange(0,n_tr,1):#What de-excited states are available to it?
				ch = F_n[j,0]#get the lower F state of the transition
				if ch == F_g:#find all the transitions with F_g as the upper
					al.append(F_n[j,:])#append properties of transition
					n_trr = n_trr + 1
			al = n.reshape(al,(n_trr,8))
			al[:,7] = al[:,7]/n.sum(al[:,7])#normalize transition strengths
			r_num = r.uniform(0,1) #generate a number between 0 and 1
			tot = al[0,7]
			low = 0
			for jj in n.arange(0,n_trr,1):
				if low<=r_num<=tot:
					F_g = al[jj,1]
					photon = al[jj,4] + r.uniform(-al[jj,5],al[jj,5])
				else:
					low = tot
					if jj == n_trr-1:
						tot = 1.0
					else:
						tot = tot+al[jj+1,7]
			excc = 0 #deexcite atom state
			ph_d = ph_d + 1.0
			photons.append(photon) #store the photon
			#print pos/(d+d_LCR)
			#-------------atom is now deexcited-------------------------------			
			#this while-loop will go on until the atom leaves the LCR
	count_t = len(photons)#total number of photons so far
	cts_dac.append(ph_d)
	cts.append(count_t)

cts_dac=n.array(cts_dac)
dacs = n.array(dacs)
nbins = int(l)

with open('Output.txt','w') as txt_file:
	for i in cts_dac:
		txt_file.write('%f\n' % (i))

with open('Params.txt','w') as txt_file:
	for param in vec:
		txt_file.write('%f\n' % (param))
#spec = n.histogram(photons,bins=nbins)
#plt.plot(spec[1][0:nbins],spec[0])

	
