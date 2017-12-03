offset = 23938.89578#laser frequency in cm^-1
lam = (1.0/offset)/100 #wl in m
En_offset = 2.0*n.pi*hbar*c/lam#offset energy in eV

data = n.loadtxt(open('FILENAME',"r"))#import the datafile with name FILENAME
exp_x = data[:,0]
exp_y = data[:,1]

exp_x_conv = exp_x*10**6#convert to Hz
exp_x_conv = exp_x_conv*2.0*n.pi*hbar#convert to eV
exp_x_conv = exp_conv+En_offset #recenter the whole thing

index = 0
for ch in exp_y:#find the highest peak in the experimental data
	if ch == n.max(exp_y):
		max_in_exp = index
		break
	index = index+1.0

index = 0
for ch in sim_y:#find the highest peak in the experimental data
	if ch == n.max(sim_y):
		max_in_sim = index
		break
	index = index+1.0


diff = exp_x_conv[max_in_exp] - sim_x[max_in_sim]#difference in locations of highest peak

exp_x_conv = exp_x_conv -diff