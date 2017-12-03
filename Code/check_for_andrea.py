ch = 0
plt.clf()
powers = []
powers_chi = []
init_powers = [0.1,0.2,0.3]
for power_plot in init_powers:
	powers.append(power_plot)
	label_plot = str(power_plot) + ' mW'
	vec[10] = 25
	vec[12] = power_plot
	vec[14] = 3
	execfile('OPV4.py')
	plt.plot(xdata,spec_sum,label=label_plot)

plt.legend()
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel('Energy (eV)')
plt.ylabel('Relative Photon Count')
plt.savefig('Power_comparison(0.1-0.3mW)_ca.png')
plt.clf()

powers_chi = []
init_powers = [0.3]
for power_plot in init_powers:
	powers.append(power_plot)
	label_plot = str(power_plot) + ' mW'
	vec[10] = 25
	vec[12] = power_plot
	vec[14] = 3
	execfile('OPV4.py')
	plt.plot(xdata,spec_sum,label=label_plot)

plt.legend()
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel('Energy (eV)')
plt.ylabel('Relative Photon Count')
plt.savefig('Power_comparison(0.3mW)_ca.png')
plt.clf()

powers_chi = []
init_powers = [0.3,1,2,3]
for power_plot in init_powers:
	powers.append(power_plot)
	label_plot = str(power_plot) + ' mW'
	vec[10] = 25
	vec[12] = power_plot
	vec[14] = 3
	execfile('OPV4.py')
	plt.plot(xdata,spec_sum,label=label_plot)

plt.legend()
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel('Energy (eV)')
plt.ylabel('Relative Photon Count')
plt.savefig('Power_comparison(a0.3-3mW)_ca.png')
plt.clf()