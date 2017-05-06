from Tkinter import *

import os.path
from math import sqrt
from math import pi
from numpy import *

def factorial(n):    if n == 1 or n==0:        return 1    else:        return n * factorial(n-1)

def binom(n, m):    b = [0] * (n + 1)    b[0] = 1    for i in xrange(1, n + 1):        b[i] = 1        j = i - 1        while j > 0:            b[j] += b[j - 1]            j -= 1    return b[m]

def angdelta(a,b,c):
    return sqrt(float(factorial((a+b-c)/2)*factorial((a-b+c)/2)*factorial((-a+b+c)/2))/float(factorial((a+b+c)/2+1)))

def sixj(a,b,c,d,e,f):
	if abs(a-b)>c or abs(c-d)>e or abs(a-e)>f or abs(b-d)>f or a+b<c or c+d<e or a+e<f or b+d<f: return 0
	outfactors = angdelta(a,e,f)*angdelta(b,d,f)*angdelta(c,d,e)/angdelta(a,b,c)
	nlo = max([(a+b+c)/2,(c+d+e)/2,(b+d+f)/2,(a+e+f)/2])
	nhi = min([(a+b+d+e)/2,(b+c+e+f)/2,(a+c+d+f)/2])
	sum = 0.0
	for n in range(nlo,nhi+1):
	    sumterm = ((-1)**n)*binom(n+1,n-(a+b+c)/2)*binom((a+b-c)/2,n-(c+d+e)/2)*binom((a-b+c)/2,n-(b+d+f)/2)*binom((b-a+c)/2,n-(a+e+f)/2)
	    sum=sum+sumterm
	return sum * outfactors

def intensities(fl,fu,nucspin,jl,ju):
	sum0=0.0; sum2=0.0; gf0=zeros(len(fl),float); fg=zeros(len(fl),float); gf2=zeros(len(fl),float); intensities=zeros([len(fl),3],float)
	for peak in range(len(fl)):
		flower=fl[peak]
		fupper=fu[peak]
		w2=(sixj(int(2*ju),int(2*fupper),int(2*nucspin),int(2*flower),int(2*jl),int(2)))**2
		gf0[peak]=(2*flower+1)*(2*fupper+1)*(2*jl+1)*w2
		fg[peak]=1.0-float((2*ju+1)*gf0[peak])/float((2*fupper+1)*(2*jl+1))
		wf=sixj(int(2*flower),int(2*fupper),2,4,2,int(2*fupper))
		wg=sixj(int(2*nucspin),int(2*fupper),int(2*ju),4,int(2*ju),int(2*fupper))
		wj=sixj(int(2*jl),int(2*ju),2,4,2,int(2*ju))*(2*flower+1)*(2*fupper+1)*(2*fupper+1)*(2*jl+1)*(2*ju+1)
		phase=(-1.0)**(nucspin-flower-jl)
		gf2[peak]=1.5*phase*wf*wg*wj*w2
		sum0 += gf0[peak]
		sum2 += gf0[peak] + gf2[peak]
	for peak in range(len(fl)):
		gf2[peak]=(gf0[peak]+gf2[peak])*sum0/sum2
		intensities[peak][0]=gf0[peak]/sum0
		intensities[peak][1]=gf2[peak]/sum0
		intensities[peak][2]=fg[peak]
	return intensities

def alpha(f,i,j):
	if i==0 or j==0: return 0
	else: return float(f*(f+1) - i*(i+1) - j*(j+1))/2.0

def beta(f,i,j):
	k = f*(f+1) - i*(i+1) - j*(j+1)
	if i==0 or i==0.5 or j==0 or j==0.5: return 0
	else: return float(3*k*(k+1) - 4*i*(i+1)*j*(j+1))/float(8*i*(2*i-1)*j*(2*j-1))

def transitions(jl,ju,ns):
	lowerfs=[]; lowerf=abs(jl-ns)
	while lowerf <= jl+ns :
		lowerfs.append(lowerf)
		lowerf = lowerf + 1
	upperfs=[]; upperf=abs(ju-ns)
	while upperf <= ju+ns :
		upperfs.append(upperf)
		upperf = upperf + 1
	transitions=[]
	for fl in lowerfs:
		if fl-1 in upperfs: transitions.append([fl,fl-1,alpha(fl,ns,jl),beta(fl,ns,jl),alpha(fl-1,ns,ju),beta(fl-1,ns,ju),0,0,0])
		if fl in upperfs and fl != 0: transitions.append([fl,fl,alpha(fl,ns,jl),beta(fl,ns,jl),alpha(fl,ns,ju),beta(fl,ns,ju),0,0,0])
		if fl+1 in upperfs: transitions.append([fl,fl+1,alpha(fl,ns,jl),beta(fl,ns,jl),alpha(fl+1,ns,ju),beta(fl+1,ns,ju),0,0,0])
	fl=[row[0] for row in transitions]
	fu=[row[1] for row in transitions]
	for peak in range(len(fl)):
		transitions[peak][6]=intensities(fl,fu,ns,jl,ju)[peak][0]
		transitions[peak][7]=intensities(fl,fu,ns,jl,ju)[peak][1]
		transitions[peak][8]=intensities(fl,fu,ns,jl,ju)[peak][2]
	return transitions

def transfreqs(t,au,bu,al,bl,w):
	transfreqs=[]
	for peak in range(shape(t)[0]):
		frequency = w + au*t[peak][4] + bu*t[peak][5] - al*t[peak][2] - bl*t[peak][3]
		transfreqs.append([t[peak][0],t[peak][1], t[peak][6], t[peak][7], t[peak][8], frequency])
	return transfreqs
		
def spectrum(x,fwhm,intensity,w,i):
	y=0.0
	for peak in range(len(w)):
		y=y+float(i[peak])*intensity*(fwhm**2/(4*((x-float(w[peak]))**2+(fwhm/2.0)**2)))
	return y

def slashform(dec):
	if int(2*dec) % 2 == 1: slashform=str(int(dec*2))+'/2'
	else: slashform=str(int(dec))
	return slashform
	
####end of racah, frequency and transition functions

jl=0;ju=0
while ((ju-jl)!=0 and (ju-jl)!=1 and (jl-ju)!=1) or (jl==0 and ju==0):
	jl = float(raw_input("jl?"))
	ju = float(raw_input("ju?"))
fwhm = int(raw_input("FWHM?"))
hyprange = int(raw_input("Maximum value of hyperfine coefficients (MHz)?"))
loadedsetfilename=raw_input("Data set file name?")
if loadedsetfilename=='' or os.path.isfile(loadedsetfilename)==False: datasetloaded=False
else:
	datasetloaded=True
	loadedset = open(loadedsetfilename,"r")
	loadedsetlines = loadedset.readlines()
	loadedset.close()
	loadedsetpoints=len(loadedsetlines)-loadedsetlines.count('\n')
	lsx=zeros(loadedsetpoints)
	lsy=zeros(loadedsetpoints)
	n=0
	for line in loadedsetlines:
		if line!="\n":
			xdata, ydata = line.split()
			lsx[n]=float(xdata)
			lsy[n]=float(ydata)
			n=n+1
	lsy=1+lsy-min(lsy)
	lsy=lsy/max(lsy)

# alter configuration below +++++++++++

canvaswidth=800.0
canvasheight=400.0
magno=3
step=int(int(fwhm)/10.0)
if step == 0: step = 1

# alter configuration above ++++++++++++
	
root=Tk()
root.title("Spectrum Viewer")

def plot(event):
	global locking
	locking=False
	global t
	flmax = max([row[0] for row in t])
	if flmax != jl+spinval.get() : t=transitions(jl,ju,spinval.get())
	w.delete("calclines")
	scale=scaleslider.get()
	mhzperpixel=scale/(canvaswidth/2)
	if ju>0: au=auval.get()
	else: au=0
	if ju>0.5: bu=buval.get()
	else: bu=0
	if jl>0: al=alval.get()
	else:al=0
	if jl>0.5: bl=blval.get()
	else:bl=0
	global f
	f=transfreqs(t,au,bu,al,bl,0)
	if fwhm*2*magno*shape(f)[0] < max([abs(row[5]) for row in f])+2*magno*fwhm:
		#if most of the spectrum is space
		nonzeroset=[]
		for peak in range(shape(f)[0]):
			if -scale-magno*fwhm<f[peak][5] and f[peak][5]<scale+magno*fwhm:
				for n in range(-magno*fwhm,magno*fwhm,step):
					m=f[peak][5]+n
					nonzeroset.append(m)
		if nonzeroset != []:
			nonzerosetsorted=sort(nonzeroset)
			a=nonzerosetsorted[0]/mhzperpixel+canvaswidth/2
			b=canvasheight-spectrum(nonzerosetsorted[0],fwhm,canvasheight,[row[5] for row in f],[row[3] for row in f])
			for x in nonzerosetsorted:
				c=x/mhzperpixel+canvaswidth/2
				d=canvasheight-spectrum(x,fwhm,canvasheight,[row[5] for row in f],[row[3] for row in f])
				id = w.create_line(a,b,c,d,tags="calclines")
				a=c;b=d
	else:
		#if the peaks are squashed up
		if scale<-(min([row[5] for row in f])+magno*fwhm): m=scale
		else: m=-min([row[5] for row in f])+magno*fwhm
		if scale<(max([row[5] for row in f])+magno*fwhm): n=scale
		else: n=max([row[5] for row in f])+magno*fwhm
		a=-m/mhzperpixel+canvaswidth/2
		b=canvasheight-spectrum(-m,fwhm,canvasheight,[row[5] for row in f],[row[3] for row in f])
		for x in range(int(-m),int(n),step):
			c=x/mhzperpixel+canvaswidth/2
			d=canvasheight-spectrum(x,fwhm,canvasheight,[row[5] for row in f],[row[3] for row in f])
			id = w.create_line(a,b,c,d,tags="calclines")
			a=c;b=d
def plotdata(event):
	global locking
	locking=False
	if datasetloaded==True:
		w.delete("datasetlines")
		scale=scaleslider.get()
		mhzperpixel=scale/(canvaswidth/2)
		a=(lsx[0]-centroidslider.get())/mhzperpixel+canvaswidth/2;b=d=canvasheight*(1-lsy[0])
		for n in range(1,len(lsx)):
			c=(lsx[n]-centroidslider.get())/mhzperpixel+canvaswidth/2
			d=canvasheight*(1-lsy[n])
			if abs(c-a)*mhzperpixel<100: id = w.create_line(a,b,c,d,fill='cyan',tags="datasetlines")
			a=c;b=d
			
def reqfreq(event):
	global locking
	global lockmatrix
	global lockfreqs
	global ausoln
	global alsoln
	global busoln
	global blsoln
	scale=scaleslider.get()
	mhzperpixel=scale/(canvaswidth/2)
	dim = 1-scalea.get()-scaleb.get()
	if jl>0 and spinval.get()>0:
		dim = dim + 1
		if jl>0.5 and spinval.get()>0.5: dim = dim +1
	if ju>0 and spinval.get()>0:
		dim = dim + 1
		if ju>0.5 and spinval.get()>0.5: dim = dim +1
	if len(w.gettags("redpeak"))==0:
		locking=False
	else:
		if locking==False:
			locking=True
			lockmatrix=zeros([dim,dim])
			lockfreqs=zeros([dim])
			ausoln=False
			alsoln=False
			busoln=False
			blsoln=False
		peak=int(w.gettags("redpeak")[3])
		w.itemconfigure("redpeak",fill="green")
		for n in xrange(dim):
			if lockmatrix[n,0]==0:
				lockfreqs[n]=float((event.x-canvaswidth/2)*mhzperpixel)
				m=0
				lockmatrix[n,m]=1.0
				m=m+1
				if scalea.get()==1 and spinval.get()>0:
					lockmatrix[n,m] = t[peak][4]-ascalefactor.get()*t[peak][2]
					ausoln=True
					m=m+1
				else:
					if ju>0 and spinval.get()>0:					
						lockmatrix[n,m] = t[peak][4]
						ausoln=True
						m=m+1
					if jl>0 and spinval.get()>0:
						lockmatrix[n,m] = -t[peak][2]
						alsoln=True
						m=m+1
				if scaleb.get()==1 and spinval.get()>0:
					lockmatrix[n,m] = t[peak][5]-bscalefactor.get()*t[peak][3]
					busoln=True
				else:
					if ju>0.5 and spinval.get()>0.5:
						lockmatrix[n,m] = t[peak][5]
						busoln=True
						m=m+1
					if jl>0.5 and spinval.get()>0.5:
						lockmatrix[n,m] = -t[peak][3]
						blsoln=True
				break
		if n==dim-1:
			if linalg.det(lockmatrix)==0:
				print "Singular matrix - repeated interval?"
			else:
				solution=linalg.solve(lockmatrix,lockfreqs)
				centroidslider.set(centroidslider.get()+solution[0])
				solnumber=1
				if ausoln==True:
					auslider.set(solution[solnumber])
					solnumber=solnumber+1
				if alsoln==True:
					alslider.set(solution[solnumber])
					solnumber=solnumber+1
				if busoln==True:
					buslider.set(solution[solnumber])
					solnumber=solnumber+1
				if blsoln==True:
					blslider.set(solution[solnumber])
			locking=False

def showsticks(event):
	w.delete("sticks")
	scale=scaleslider.get()
	mhzperpixel=scale/(canvaswidth/2)
	for peak in range(shape(f)[0]):
		xpos=f[peak][5]/mhzperpixel+canvaswidth/2
		if abs(event.x-xpos)<5 and abs(event.y-(canvasheight*(1-f[peak][3])))<5:
			w.create_line(xpos,canvasheight,xpos,canvasheight*(1-f[peak][3]), arrow=LAST, fill='red',tags=("calclines","sticks","redpeak",peak))
			peakinfo=slashform(f[peak][0])+"->"+slashform(f[peak][1])+"\n frequency = "+str(f[peak][5])+"\n intensity@90 = "+str(f[peak][3])
			w.create_text(xpos,canvasheight-400*f[peak][3]-20,text=peakinfo,tags=("calclines","sticks"))
		else:
			w.create_line(xpos,canvasheight,xpos,canvasheight*(1-f[peak][3]), arrow=LAST, tags=("calclines","sticks"))

w=Canvas(root,background="white",width=canvaswidth,height=canvasheight)
w.pack(side='top')
w.bind('<Shift-Button-1>',reqfreq)
w.bind('<Button-1>',showsticks)

v=Canvas(root,background="white",width=canvaswidth,height=30)
v.pack(side='top')

scaleframe=Frame(root)
scaleframe.pack(side='top')
def tickmarksandlabels(scale):
	id = v.create_line(0,5,canvaswidth,5)
	mhzperpixel=scale/(canvaswidth/2)
	spacings=[10,20,50,100,200,500,1000,2000,5000,10000,20000]
	for spacing in spacings:
		if scale/spacing < 10: break
	for tickmark in range(scale/spacing):
		id=v.create_line(canvaswidth/2+tickmark*spacing/mhzperpixel,5,canvaswidth/2+tickmark*spacing/mhzperpixel,15)
		id=v.create_text(canvaswidth/2+tickmark*spacing/mhzperpixel,25,text=tickmark*int(spacing))
		id=v.create_line(canvaswidth/2-tickmark*spacing/mhzperpixel,5,canvaswidth/2-tickmark*spacing/mhzperpixel,15)
		id=v.create_text(canvaswidth/2-tickmark*spacing/mhzperpixel,25,text=-tickmark*int(spacing))
def rescale(event):
	v.delete(ALL)
	scale=scaleslider.get()
	tickmarksandlabels(scale)
	plotdata(event)
	plot(event)
scaletext=Label(scaleframe,text='horizontal scale')
scaletext.pack(side='left')
scaleslider=IntVar()
scaleslider.set(100)
scaleslider=Scale(scaleframe,orient='horizontal',from_=100,to=50000, showvalue=0,length=600,width=12,resolution=100,command=rescale)
scaleslider.pack(side='left')

optionsframe=Frame(root)
optionsframe.pack(side='top')

def spinchange(event):
	spintemp=spinval.get()
	spin=int(2*spintemp)/2.0
	spinval.set(spin)
	global t
	t=transitions(jl,ju,spin)
	plot(event)
spinframe=Frame(optionsframe)
spinframe.pack(side='left')
spintext=Label(spinframe, text='I=')
spintext.pack(side='left')
spinval=DoubleVar()
spinval.set(1.0)
spinval_entry=Entry(spinframe, width=4, relief=RIDGE, textvariable=spinval)
spinval_entry.pack(side='left')
spinval_entry.bind('<Return>', spinchange)
spinval_entry.bind('<FocusOut>', spinchange)

scalea=IntVar()
scalea.set(0)
if jl>0 and ju>0:
	def scalinga(event):
		if scalea.get()==1:
			alval.set(auval.get()*ascalefactor.get())
			alslider.set(auval.get()*ascalefactor.get())
		else: plot(event)
	scaleaframe=Frame(optionsframe)
	scaleaframe.pack(side='left',padx=50)
	scaleabutton=Checkbutton(scaleaframe, text='Scale Al/Au=',variable=scalea)
	scaleabutton.pack(side='left')
	ascalefactor=DoubleVar()
	ascalefactor.set('1.0')
	ascalefactor_entry=Entry(scaleaframe, width=6, relief=RIDGE, textvariable=ascalefactor)
	ascalefactor_entry.pack(side='left')
	ascalefactor_entry.bind('<Return>', scalinga)
	ascalefactor_entry.bind('<FocusOut>', scalinga)
	def s(event):
		alval.set(auval.get()*ascalefactor.get())
		alslider.set(auval.get()*ascalefactor.get())
	scaleabutton.bind('<1>', s)

scaleb=IntVar()
scaleb.set(0)
if jl>0.5 and ju>0.5:
	def scalingb(event):
		if scaleb.get()==1:
			blval.set(buval.get()*bscalefactor.get())
			blslider.set(buval.get()*bscalefactor.get())
		else: plot(event)
	scalebframe=Frame(optionsframe)
	scalebframe.pack(side='left',padx=50)
	scalebbutton=Checkbutton(scalebframe, text='Scale Bl/Bu=',variable=scaleb)
	scalebbutton.pack(side='left')
	bscalefactor=DoubleVar()
	bscalefactor.set('1.0')
	bscalefactor_entry=Entry(scalebframe, width=6, relief=RIDGE, textvariable=bscalefactor)
	bscalefactor_entry.pack(side='left')
	bscalefactor_entry.bind('<Return>', scalingb)
	bscalefactor_entry.bind('<FocusOut>', scalingb)
	def t(event):
		blval.set(buval.get()*bscalefactor.get())
		blslider.set(buval.get()*bscalefactor.get())
	scalebbutton.bind('<1>', t)

def dump():
	spin=spinval.get()
	print "jl = %.1f , ju = %.1f , i = %.1f" % (jl,ju,spin)
	if ju>0 and spin>0: 
		au=auval.get()
		print "Au = %f" % (au)
	else: au=0
	if ju>0.5 and spin>0.5:
		bu=buval.get()
		print "Bu = %f" % (bu)
	else: bu=0
	if jl>0 and spin>0: 
		al=alval.get()
		print "Al = %f" % (al)
	else: al=0
	if jl>0.5 and spin>0.5: 
		bl=blval.get()
		print "Bl = %f" % (bl)
	else: bl=0
	t=transitions(jl,ju,spin)
	f=transfreqs(t,au,bu,al,bl,0)
	print "FWHM = "+str(fwhm)+"\n" 
 	print"fl   fu    alpha_l      beta_l     alpha_u     beta_u       freq      int@90  pumpfrac"
	print"==   ==    =======      ======     =======     ======       ========  ======  ========"
	for r in range(shape(t)[0]):
		print " %-3s  %-3s %10.6f  %10.6f  %10.6f  %10.6f  %10.2f  %.4f  %.4f" % (slashform(t[r][0]),slashform(t[r][1]),t[r][2],t[r][3],t[r][4],t[r][5],f[r][5],t[r][7],t[r][8])
	savename=raw_input("Save spectrum as?")
	if savename!='': 
		if os.path.isfile(savename)==True: print "File already exists"
		else: 
			outputfile = open(savename,"w")
			for freq in range(int(min([row[5] for row in f])-magno*fwhm),int(max([row[5] for row in f])+magno*fwhm),step):
				if spectrum(freq,fwhm,canvasheight,[row[5] for row in f],[row[3] for row in f])>0.01:
					outputline=str(freq)+"    "+str(spectrum(freq,fwhm,100,[row[5] for row in f],[row[3] for row in f]))+"\n"
					outputfile.write(outputline)
			outputfile.close()
dumpbutton=Button(optionsframe, text='Dump', command=dump)
dumpbutton.pack(side='left',padx=50)

t=transitions(jl,ju,spinval.get())
locking=False

if datasetloaded==True:
	centroidframe=Frame(root)
	centroidframe.pack(side='top')
	def centroidslide(event):
		centroidval.set(centroidslider.get())
		plotdata(event)
	def centroidtyped(event):
		centroidslider.set(centroidval.get())
		plotdata(event)
	centroidtext=Label(centroidframe,text=' w ')
	centroidtext.pack(side='left')
	centroidslider=Scale(centroidframe,orient='horizontal',from_=min(lsx),to=max(lsx),length=600,width=12,showvalue=0,resolution=0.01,command=centroidslide)
	centroidslider.pack(side='left')
	centroidval=DoubleVar()
	centroidval.set((min(lsx)+max(lsx))/2.0)
	centroid_entry=Entry(centroidframe, width=6, relief=RIDGE, textvariable=centroidval)
	centroid_entry.pack(side='left')
	centroid_entry.bind('<Return>', centroidtyped)
	centroid_entry.bind('<FocusOut>', centroidtyped)
		
if ju>0:
	auframe=Frame(root)
	auframe.pack(side='top')
	def auslide(event):
		auval.set(auslider.get())
		if jl>0: scalinga(event)
		else: plot(event)
	def autyped(event):
		auslider.set(auval.get())
		if jl>0: scalinga(event)
		else: plot(event)
	autext=Label(auframe, text='Au ')
	autext.pack(side='left')
	auslider=Scale(auframe,orient='horizontal',from_=-hyprange,to=hyprange,length=600,width=12,showvalue=0,resolution=0.01,command=auslide)
	auslider.pack(side='left')
	auval=DoubleVar()
	auval.set('0.0')
	auval_entry=Entry(auframe, width=6, relief=RIDGE, textvariable=auval)
	auval_entry.pack(side='left')
	auval_entry.bind('<Return>', autyped)
	auval_entry.bind('<FocusOut>', autyped)
	
if ju>0.5:
	buframe=Frame(root)
	buframe.pack(side='top')
	def buslide(event):
		buval.set(buslider.get())
		if jl>0.5: scalingb(event)
		else: plot(event)
	def butyped(event):
		buslider.set(buval.get())
		if jl>0.5: scalingb(event)
		else: plot(event)
	butext=Label(buframe, text='Bu ')
	butext.pack(side='left')
	buslider=Scale(buframe,orient='horizontal',from_=-hyprange,to=hyprange,length=600,width=12,showvalue=0,resolution=0.01,command=buslide)
	buslider.pack(side='left')
	buval=DoubleVar()
	buval.set('0.0')
	buval_entry=Entry(buframe, width=6, relief=RIDGE, textvariable=buval)
	buval_entry.pack(side='left')
	buval_entry.bind('<Return>', butyped)
	buval_entry.bind('<FocusOut>', butyped)
	
if jl>0:
	alframe=Frame(root)
	alframe.pack(side='top')
	def alslide(event):
		alval.set(alslider.get())
		if ju>0: scalinga(event)
		plot(event)
	def altyped(event):
		alslider.set(alval.get())
		if ju>0: scalinga(event)
		plot(event)
	altext=Label(alframe, text='Al ')
	altext.pack(side='left')
	alslider=Scale(alframe,orient='horizontal',from_=-hyprange,to=hyprange,length=600,width=12,showvalue=0,resolution=0.01,command=alslide)
	alslider.pack(side='left')
	alval=DoubleVar()
	alval.set('0.0')
	alval_entry=Entry(alframe, width=6, relief=RIDGE, textvariable=alval)
	alval_entry.pack(side='left')
	alval_entry.bind('<Return>', altyped)
	alval_entry.bind('<FocusOut>', altyped)
	
if jl>0.5:
	blframe=Frame(root)
	blframe.pack(side='top')
	def blslide(event):
		blval.set(blslider.get())
		if ju>0.5: scalingb(event)
		plot(event)
	def bltyped(event):
		blslider.set(blval.get())
		if ju>0.5: scalingb(event)
		plot(event)
	bltext=Label(blframe, text='Bl ')
	bltext.pack(side='left')
	blslider=Scale(blframe,orient='horizontal',from_=-hyprange,to=hyprange,length=600,width=12,showvalue=0,resolution=0.01,command=blslide)
	blslider.pack(side='left')
	blval=DoubleVar()
	blval.set('0.0')
	blval_entry=Entry(blframe, width=6, relief=RIDGE, textvariable=blval)
	blval_entry.pack(side='left')
	blval_entry.bind('<Return>', bltyped)
	blval_entry.bind('<FocusOut>', bltyped)

root.mainloop()
