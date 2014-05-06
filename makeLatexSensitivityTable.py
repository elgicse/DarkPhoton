import os
pwd = os.path.dirname(os.path.abspath(__file__))

lol = []

def zip_lists(*args):
	return map(list, zip(*args))


def latex_float(fl):
	if fl == 0. or fl == -0.:
		return "0.0"
	float_str = "{0:.1e}".format(fl)
	if "e" in float_str:
		base, exponent = float_str.split("e")
		if float(exponent) == 0.:
			return r"{0}".format(float(base))
		elif float(exponent) == 1.:
			return "%s"%(int(fl))
		else:
			return r"${0} \times 10^{{{1}}}$".format(float(base), int(exponent))
	else:
		return float_str


for root, dirs, filenames in os.walk(pwd):
    for f in filenames:
    	if "sensitivities" in f and "~" not in f:
    		#print f
    		l = []
    		l.append(os.path.join(root, f))
    		f1 = f.replace("decays","_").replace(".txt","").replace("-","")
    		for t in f1.split("_"):
    			try:
    				l.append(float(t))
    			except ValueError:
    				pass
    		if len(l) == 2:
    			l.extend([0.0])
    		l[1] = int(l[1])
    		l[2] = int(l[2])
    		#print l
    		lol.append(l)
zipped = zip_lists(*lol)
#print zipped
sorted_files = [list(x) for x in zip(*sorted(zipped, key=lambda pair: pair[1]))]
temp = sorted_files[:]
sorted_files = sorted(temp, key=lambda pair: pair[0])
#print sorted_files

eps, mass, prodrate, BRee, BRmumu, accee, accmumu, totalacc, nevents, ngenerated = [], [], [], [], [], [], [], [], [], []

for datafile in sorted_files:
	inputfile = open(datafile[2],'r')
	for line in inputfile:
		line = line.split()
		try:
			veps = float(line[0])
			vmass = float(line[1])
			vprodrate = float(line[2])
			vBRee = float(line[3])
			vBRmumu = float(line[4])
			vaccee = float(line[5])
			vaccmumu = float(line[6])
			vtotalacc = float(line[7])
			vnevents = float(line[8])

			eps.append(veps)
			mass.append(vmass)
			prodrate.append(vprodrate)
			BRee.append(vBRee)
			BRmumu.append(vBRmumu)
			accee.append(vaccee)
			accmumu.append(vaccmumu)
			totalacc.append(vtotalacc)
			nevents.append(int(vnevents))
			ngenerated.append(int(datafile[1]))
		except:
			continue

#for brmu in BRmumu:
#	if brmu < 0:
#		brmu = -1*brmu

alldata = sorted(zip(mass, eps, prodrate, BRee, BRmumu, accee, accmumu, totalacc, ngenerated, nevents), key=lambda x: (x[0], x[1]))

latexhdr = r"\begin{tabular}{llllllllll}\hline $M_{\gamma'}$ [GeV] & $\varepsilon$ & Prod. rate ($\chi$) & $BR(\gamma'\rightarrow ee)$ & $BR(\gamma'\rightarrow \mu\mu)$ & Acc.($ee$) & Acc.($\mu\mu$) & Total acc. & Gen. sample & Events in SHiP\\\hline"
print latexhdr

for dataline in alldata:
	if dataline[9] > 0:
		latexline = []
		for datum in dataline:
			latexline.append(latex_float(datum))
		separator = r" & "
		print separator.join(latexline) + r"\\"

latexftr = r"\hline\end{tabular}"
print latexftr