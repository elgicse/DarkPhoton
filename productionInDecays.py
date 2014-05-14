# Check $ROOTSYS/tutorials/pythia/pythia8.C
# for example program
import ROOT as r

def loadLibs():
	r.gSystem.Load("$PYTHIA8/lib/libpythia8")
	r.gSystem.Load("libEG")
	r.gSystem.Load("libEGPythia8")

if __name__ == '__main__':

	loadLibs()
	pythia = r.TPythia8()
	pythia.ReadString("SoftQCD:nonDiffractive = on")
	pythia.ReadString("Beams:idA = 2212")
	pythia.ReadString("Beams:idB = 2212")
	pythia.ReadString("Beams:frameType = 2")
	pythia.ReadString("Beams:eA = 400.")
	pythia.ReadString("Beams:eB = 0.")
	pyt = pythia.Pythia8()
	pyt.init()
	nev = 10
	for i in xrange(nev):
		pythia.GenerateEvent()
		seen = []
		print "Event %s: "%i, pyt.info.name()
		names = [event[part].name() for part in pyt.event[0].daughterList()]
		if i ==4:
			print " ".join(childrenNames[:])
		for p in xrange(pyt.event.size()):
		    part = pyt.event[p]
		    if (part.id() == 22) and (p not in seen):
		    	mum = pyt.event[pyt.event[p].mother1()]
		    	children = mum.daughterList()
		    	seen.extend(list(children))
		    	childrenList = []
		    	for c in xrange(len(children)):
		    		childrenList.append(pyt.event[children[c]])
		    	childrenNames = [child.name() for child in childrenList]
		    	print mum.name(), " ".join(childrenNames[:])
