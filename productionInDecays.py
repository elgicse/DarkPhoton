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
	pythia.ReadString("SoftQCD:singleDiffractive = on")
	pythia.ReadString("SoftQCD:doubleDiffractive = on")
	pythia.Initialize(2212, 2212, 14000.)
	nev = 10
	parts = r.TClonesArray("TParticle", 1000)

	for i in xrange(nev):
	    pythia.GenerateEvent()
	    pythia.ImportParticles(parts,"All")
	    for p in xrange(parts.GetEntriesFast()):
	        part = parts.At(p)
	        print part.GetPdgCode()
