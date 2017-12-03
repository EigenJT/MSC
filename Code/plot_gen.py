lik_absorption = []
ch = 0
plt.clf()
powers = []
powers_chi = []
power_first = n.arange(0.1,40.05,0.05)
init_powers = [0.1,1,5,10,20,30,40]
for power_plot in init_powers:
	powers.append(power_plot)
	label_plot = str(power_plot) + ' mW'
	vec[10] = 300
	vec[12] = power_plot
	vec[14] = 0.4
	execfile('OPV4.py')
	plt.plot(xdata,spec_sum,label=label_plot)
	powers_chi.append(chi_2_red)
	lik_absorption.append(F_n[:,12])

plt.legend()
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel('Energy (eV)')
plt.ylabel('Relative Photon Count')
plt.savefig('Power_comparison(1-40).png')
plt.clf()


power_second = n.arange(40.1,50.05,0.05)
sec_power = [50,60,75,90,100]
for power_plot in sec_power:
	powers.append(power_plot)
	label_plot = str(power_plot) + ' mW'
	vec[10] = 300
	vec[12] = power_plot
	vec[14] = 0.4
	execfile('OPV4.py')
	plt.plot(xdata,spec_sum,label=label_plot)
	powers_chi.append(chi_2_red)
	lik_absorption.append(F_n[:,12])

plt.legend()
plt.xlabel('Energy (eV)')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ylabel('Relative Photon Count')
plt.savefig('Power_comparison(50-100).png')
plt.clf()

plt.plot(powers,powers_chi)
plt.plot([1],powers_chi[6],'k|')
plt.xlabel('Power (mW)')
plt.ylabel('Reduced $\chi^{2}$')
plt.savefig('power_chi.png')
plt.clf()



for power_plot in [110,120,130,140,150]:
	powers.append(power_plot)
	label_plot = str(power_plot) + ' mW'
	vec[10] = 300
	vec[12] = power_plot
	vec[14] = 0.4
	execfile('OPV4.py')
	powers_chi.append(chi_2_red)
	plt.plot(xdata,spec_sum,label=label_plot)
	lik_absorption.append(F_n[:,12])

plt.legend()
plt.xlabel('Energy (eV)')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ylabel('Relative Photon Count')
plt.savefig('Power_comparison(110-150).png')
plt.clf()


for power_plot in [200,250,300,350,400]:
	powers.append(power_plot)
	label_plot = str(power_plot) + ' mW'
	vec[10] = 300
	vec[12] = power_plot
	vec[14] = 0.4
	execfile('OPV4.py')
	powers_chi.append(chi_2_red)
	plt.plot(xdata,spec_sum,label=label_plot)
	lik_absorption.append(F_n[:,12])

plt.legend()
plt.xlabel('Energy (eV)')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ylabel('Relative Photon Count')
plt.savefig('Power_comparison(200-400).png')
plt.clf()

for power_plot in [450,500,650,700,850]:
	powers.append(power_plot)
	label_plot = str(power_plot) + ' mW'
	vec[10] = 300
	vec[12] = power_plot
	vec[14] = 0.4
	execfile('OPV4.py')
	powers_chi.append(chi_2_red)
	plt.plot(xdata,spec_sum,label=label_plot)
	lik_absorption.append(F_n[:,12])

plt.legend()
plt.xlabel('Energy (eV)')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ylabel('Relative Photon Count')
plt.savefig('Power_comparison(450-850).png')
plt.clf()

for power_plot in [900,950,1000,1050]:
	powers.append(power_plot)
	label_plot = str(power_plot) + ' mW'
	vec[10] = 300
	vec[12] = power_plot
	vec[14] = 0.4
	execfile('OPV4.py')
	powers_chi.append(chi_2_red)
	plt.plot(xdata,spec_sum,label=label_plot)
	lik_absorption.append(F_n[:,12])

plt.legend()
plt.xlabel('Energy (eV)')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ylabel('Relative Photon Count')
plt.savefig('Power_comparison(900-1050).png')
plt.clf()


deleted = 0
val = 0

len_abs = len(powers)
#for check in n.arange(1,len_abs,1):
#	val = lik_absorption[check-deleted]
#	valm1 = lik_absorption[check-1-deleted]
#	if n.sum(val-valm1)==0:
#		del lik_absorption[check-deleted]
#		deleted = deleted+1

lik_absorption=n.reshape(lik_absorption,(len_abs-deleted,6))
print 'Power plots done'



ch = 0
temp_chi = []
#temps = n.arange(100,2010,10)
temps = [100,150,200,300,600,1000,1500,2000]
for temp_plot in temps:
	label_plot = str(temp_plot) + 'K'
	vec[12] = 1
	vec[10] = temp_plot
	vec[14] = 0.4
	execfile('OPV4.py')
	temp_chi.append(chi_2_red)
	print temp_plot
	plt.plot(xdata,spec_sum,label=label_plot)
	ch=ch+1

plt.legend()
plt.xlabel('Energy (eV)')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ylabel('Relative Photon Count')
plt.savefig('temp_comparison.png')
plt.minorticks_on()
plt.clf()

min = n.argsort(temp_chi)[0]
plt.plot(temps,temp_chi)
plt.plot([300],temp_chi[3],'ko',label = 'T$_{est}$ $\sim$ 300 K, $\chi^2$ = %f'%(temp_chi[3]))
plt.plot(temps[min],temp_chi[min],'ro',label = 'T$_{min}$ $\sim$ %i K, $\chi^2$ = %f'%(temps[min],temp_chi[min]))
plt.xlabel('Temperature (K)')
plt.ylabel('Reduced $\chi^{2}$')
plt.legend(numpoints=1)
plt.minorticks_on()
plt.savefig('temp_chi.png')
plt.clf()
print 'Temp plot done'

dist_chi = []
dists = [0.1,0.5,0.75,1,2,5]
#dists = n.arange(0.05,5.05,0.05)
for dist_plot in dists:
	label_plot = str(dist_plot) + ' m'
	vec[12] = 1
	vec[10] = 300
	vec[14] = dist_plot
	execfile('OPV4.py')
	dist_chi.append(chi_2_red)
	plt.plot(xdata,spec_sum,label=label_plot)

plt.legend()
plt.xlabel('Energy (eV)')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ylabel('Relative Photon Count')
plt.savefig('dist_comparison.png')
plt.clf()

plt.plot(dists,dist_chi)
plt.plot([0.4],dist_chi[7],'k|')
plt.xlabel('Distance (m)')
plt.ylabel('Reduced $\chi^{2}$')
plt.savefig('dist_chi.png')
plt.clf()

print 'Dist plot done'


for abs_plot in n.arange(0,n_tr,1):
	label_upper = str(F_n[abs_plot,0])
	label_lower = str(F_n[abs_plot,1])
	label_str = 'F$_{\mathrm{upper}}$='+label_upper+'$\leftrightarrow$'+'F$_{\mathrm{lower}}$='+label_lower
	plt.semilogx(powers,lik_absorption[:,abs_plot],label=label_str)
	plt.ylabel('Absorbtion Likelihood')
	plt.xlabel('Power (mW)')

plt.legend(loc='upper left')
plt.savefig('abs.png')
print 'All done!'

