lik_absorption = []

ch = 0
for power_plot in [0.001,0.005,0.01,0.05,0.1,0.5,1]:
	label_list = ['1','10','10$^2$','10$^3$','5x10$^3$','10$^4$','5x10$^4$','10$^5$']
	label_plot = str(power_plot) + ' mW'
	vec[10] = 2000
	vec[12] = power_plot
	vec[14] = 0.4
	execfile('OPV4.py')
	print F_n[:,12]
	lik_absorption.append(F_n[:,12])

plt.legend()
plt.xlabel('Energy (eV)')
plt.ylabel('Relative Photon Count')
plt.savefig('Power_comparison(0_001-1).png')
plt.clf()


print 'Power plots done'



ch = 0
for temp_plot in [2,n.log(2.5)+2*n.log(10),n.log(5)+2*n.log(10),n.log(7.5)+2*n.log(10),3,n.log(2.5)+3*n.log(10),n.log(5)+3*n.log(10)]:
	label_list = ['10$^2$','2.5x10$^2$','5x10$^2$','7.5x10$^2$','10$^3$','2.5x10$^3$','5x10$^3$']
	label_plot = label_list[ch] + 'K'
	vec[12] = 1
	vec[10] = 10**(temp_plot)
	vec[14] = 0.4
	execfile('OPV4.py')
	ch=ch+1

plt.legend()
plt.xlabel('Energy (eV)')
plt.ylabel('Relative Photon Count')
plt.savefig('temp_comparison.png')
plt.clf()

print 'Temp plot done'

for dist_plot in [0.01,0.1,0.5,0.75,1,2,5]:
	label_plot = str(dist_plot) + 'm'
	vec[12] = 1
	vec[10] = 2000
	vec[14] = dist_plot
	execfile('OPV4.py')

plt.legend()
plt.xlabel('Energy (eV)')
plt.ylabel('Relative Photon Count')
plt.savefig('dist_comparison.png')
plt.clf()

print 'Dist plot done'

powers = [0.001,0.005,0.01,0.05,0.1,0.5,1]
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

