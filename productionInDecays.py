# Check $ROOTSYS/tutorials/pythia/pythia8.C
# for example program
from __future__ import division
import sys
import os
import ROOT as r
import numpy as np
import math
import gc

from settings import *
from functions import *
from proton_bremsstrahlung import scanPDF, computeNEvents

def loadLibs():
	r.gSystem.Load("$PYTHIA8/lib/libpythia8")
	r.gSystem.Load("libEG")
	r.gSystem.Load("libEGPythia8")
	r.gSystem.Load("$PYTHIA8/PythiaDict/pythiaDict.so")

def readSettings(pythia):
	pythia.ReadString("SoftQCD:nonDiffractive = on")
	pythia.ReadString("Beams:idA = 2212")
	pythia.ReadString("Beams:idB = 2212")
	pythia.ReadString("Beams:frameType = 2")
	pythia.ReadString("Beams:eA = 400.")
	pythia.ReadString("Beams:eB = 0.")
	pythia.ReadString("Random:setSeed = on")
	pythia.ReadString("Random:seed = 0")
	# No event record printout.
	pythia.ReadString("Next:numberShowInfo = 0")
	pythia.ReadString("Next:numberShowProcess = 0")
	pythia.ReadString("Next:numberShowEvent = 0")


def manipulatePhysics(mass, pythia):
	if (pi0Start < mass < pi0Stop):
		# use pi0 -> gamma A'
		selectedMum = 111
		pythia.ReadString("111:oneChannel = 1 1 0 22 23")
	elif etaStart < mass < etaStop:
		# use eta -> gamma A'
		selectedMum = 221
		pythia.ReadString("221:oneChannel = 1 1 0 22 23")
	elif omegaStart < mass < omegaStop:
		# use omega -> pi0 A'
		selectedMum = 223
		pythia.ReadString("223:oneChannel = 1 1 0 111 23")
	elif eta1Start < mass < eta1Stop:
		# use eta' -> gamma A'
		selectedMum = 331
		pythia.ReadString("331:oneChannel = 1 1 0 22 23")
	else:
		print "ERROR: please enter a nicer mass."
		return -1
	return selectedMum




def prepareHistograms(mass, binsp, binstheta):
	tmax = v1ThetaMax
	tmin = 0.#(-1.)*tmax # CONTROLLARE PPGAMMA
	anglestep = tmax/binstheta # CONTROLLARE PPGAMMA
	pMax = math.sqrt(protonEnergy**2.-mass**2.)
	momentumStep = pMax/binsp
	hPDF = r.TH2F("hPDF_m%s"%(mass) ,"hPDF_m%s"%(mass),
		binsp,0.5*momentumStep,pMax-0.5*momentumStep,
		binstheta,tmin-0.5*anglestep,tmax+0.5*anglestep)
	hPDF.SetTitle("PDF for A' production (m_{A'}=%s GeV)"%(mass))
	hPDF.GetXaxis().SetTitle("P_{A'} [GeV]")
	hPDF.GetYaxis().SetTitle("#theta_{A'} [rad]")
	hPDFtheta = r.TH1F("hPDFtheta_m%s"%(mass),
		"hPDFtheta_m%s"%(mass),
		binstheta,tmin-0.5*anglestep,tmax+0.5*anglestep)
	hPDFp = r.TH1F("hPDFp_m%s"%(mass),
		"hPDFp_m%s"%(mass),
		binsp,0.5*momentumStep,pMax-0.5*momentumStep)
	hPDFp.GetXaxis().SetTitle("P_{A'} [GeV]")
	hPDFtheta.GetXaxis().SetTitle("#theta_{A'} [rad]")
	return hPDF, hPDFp, hPDFtheta


def computeMom(px, py, pz):
	return math.sqrt(px*px + py*py + pz*pz)

def computeTheta(px, py, pz):
	return math.atan2(math.sqrt(px*px + py*py), pz)


def makePDF(howmany, mass, binsp=90, binstheta=40):
	gc.collect()
	loadLibs()
	pythia = r.TPythia8()
	readSettings(pythia)
	pythia.ReadString("23:m0 = %s"%mass)
	pythia.ReadString("23:name = A'")
	# Reduce width of B.W.
	pythia.ReadString("23:mWidth = 0.")
	pythia.ReadString("23:isResonance = on")
	pythia.ReadString("23:doForceWidth = on")
	# Store only the A', we will make it decay later
	pythia.ReadString("23:mayDecay = off")
	selectedMum = manipulatePhysics(mass, pythia)
	pyt = pythia.Pythia8()
	nev = howmany
	pyt.init()
	# Variables for the tree
	A_Px = np.zeros(1, dtype=float)
	A_Py = np.zeros(1, dtype=float)
	A_Pz = np.zeros(1, dtype=float)
	A_E = np.zeros(1, dtype=float)
	A_P = np.zeros(1, dtype=float)
	A_Theta = np.zeros(1, dtype=float)
	# Make or update the tree
	if os.path.isfile("out/PythiaData/mesonDecays_%s.root"%mass):
		outfile = r.TFile("out/PythiaData/mesonDecays_%s.root"%mass,"update")
		tree = outfile.Get("mesonDecays")
		tree.SetBranchAddress("A_Px", A_Px)
		tree.SetBranchAddress("A_Py", A_Py)
		tree.SetBranchAddress("A_Pz", A_Pz)
		tree.SetBranchAddress("A_E",  A_E)
		tree.SetBranchAddress("A_P",  A_P)
		tree.SetBranchAddress("A_Theta",  A_Theta)
		hPDF = outfile.Get("hPDF_m%s"%(mass))
		hPDFp = outfile.Get("hPDFp_m%s"%(mass))
		hPDFtheta = outfile.Get("hPDFtheta_m%s"%(mass))
	else:
		outfile = r.TFile("out/PythiaData/mesonDecays_%s.root"%mass,"recreate")
		tree = r.TTree("mesonDecays","mesonDecays")
		tree.Branch("A_Px", A_Px, "A_Px/D")
		tree.Branch("A_Py", A_Py, "A_Py/D")
		tree.Branch("A_Pz", A_Pz, "A_Pz/D")
		tree.Branch("A_E",  A_E , "A_E/D")
		tree.Branch("A_P",  A_P , "A_P/D")
		tree.Branch("A_Theta",  A_Theta , "A_Theta/D")
		# Histograms
		hPDF, hPDFp, hPDFtheta = prepareHistograms(mass, binsp, binstheta)
	print "Selected m(A') = %s GeV"%mass
	# Event loop
	for i in xrange(nev):
		if not (i%500):
			print "Generating event %s..."%i
			gc.collect()
		pythia.GenerateEvent()
		seen = []
		# Particle loop
		for p in xrange(pyt.event.size()):
		    part = pyt.event[p]
		    if (part.id() == 23):
		    	A_Px[0] = part.px()
		    	A_Py[0] = part.py()
		    	A_Pz[0] = part.pz()
		    	A_E[0] = part.e()
		    	mom = computeMom(part.px(), part.py(), part.pz())
		    	theta = computeTheta(part.px(), part.py(), part.pz())
		    	A_P[0] = mom
		    	A_Theta[0] = theta
		    	tree.Fill()
		    	hPDF.Fill(mom, theta)
		    	hPDFp.Fill(mom)
		    	hPDFtheta.Fill(theta)
		    	## Check event
		    	#mum = pyt.event[part.mother1()]
		    	#children = mum.daughterList()
		    	#childrenList = []
		    	#for c in xrange(len(children)):
		    	#	childrenList.append(pyt.event[children[c]])
		    	#childrenNames = [child.name() for child in childrenList]
		    	#print (mum.name() + " " + " ".join(childrenNames[:]))
		del seen
	tree.Write("",5) # TObject::kOverwrite
	hPDF.Write("",5)
	hPDFp.Write("",5)
	hPDFtheta.Write("",5)
	outfile.Close()
	gc.collect()



def makePOTs(howmany, parameters):
	""" Pass in the name of a parameters file 
	...useless for now! """
	loadLibs()
	pythia = r.TPythia8()
	readSettings(pythia)
	pyt = pythia.Pythia8()
	nev = howmany
	changedPars = parameters
	pyt.readFile(changedPars)
	pyt.init()

	# Variables for the tree
	decayName = r.TString()
	mumMass = np.zeros(1, dtype=float)
	mumId = np.zeros(1, dtype=int)
	eventNumber = np.zeros(1, dtype=int)
	childrenId = r.vector('int')()
	childrenMass = r.vector('float')()

	previousEvents = 0

	if os.path.isfile("out/PythiaData/tree.root"):
		outfile = r.TFile("out/PythiaData/tree.root","update")
		tree = outfile.Get("events")
		tree.GetEntry(tree.GetEntries()-1)
		previousEvents = int(tree.eventNumber)
		tree.SetBranchAddress("eventNumber", eventNumber)
		tree.SetBranchAddress("decayName",   decayName)
		tree.SetBranchAddress("mumId",       mumId)
		tree.SetBranchAddress("mumMass",     mumMass)
		tree.SetBranchAddress("childrenId",  childrenId)
		tree.SetBranchAddress("childrenMass",childrenMass)
	else:
		outfile = r.TFile("out/PythiaData/tree.root","recreate")
		tree = r.TTree("events","events")
		tree.Branch("eventNumber",eventNumber,"eventNumber/I")
		tree.Branch("decayName",decayName,300,0)
		tree.Branch("mumId",mumId,"mumId/I")
		tree.Branch("mumMass",mumMass,"mumMass/D")
		tree.Branch("childrenId",childrenId)
		tree.Branch("childrenMass",childrenMass)


	# Create the tree branches
	for i in xrange(nev):
		if not (i%500):
			print "Generating event %s..."%i
			gc.collect()
		pythia.GenerateEvent()
		seen = []
		names = [event[part].name() for part in pyt.event[0].daughterList()]
		for p in xrange(pyt.event.size()):
		    part = pyt.event[p]
		    if (part.id() == 22) and (part.mother1() not in seen) and (p not in seen):
		    	seen.append(part.mother1())
		    	mum = pyt.event[part.mother1()]
		    	children = mum.daughterList()
		    	# Conto un solo fotone se ce n'e' piu' d'uno nello stesso decay!
		    	seen.extend(list(children))
		    	childrenList = []
    			eventNumber[0] = previousEvents+i
		    	mumId[0] = mum.id()
		    	mumMass[0] = mum.m0()
		    	childrenId.resize(0)
		    	childrenMass.resize(0)
		    	for c in xrange(len(children)):
		    		childrenList.append(pyt.event[children[c]])
		    		childrenId.push_back(pyt.event[children[c]].id())
		    		childrenMass.push_back(pyt.event[children[c]].m0())
		    	childrenNames = [child.name() for child in childrenList]
		    	decayName.Resize(0)
		    	decayName.Append(mum.name() + " " + " ".join(childrenNames[:]))
		    	tree.Fill()

	tree.Write("",5) # TObject::kOverwrite
	outfile.Close()


# Simple analysis:
#fi = r.TFile("out/PythiaData/tree.root","read")
#tr = file.Get("events")                          
#for event in tr:
#	print event.decayName, event.mumMass, event.childrenMass[-1]
#for event in tr:
#    if event.mumMass > 0.5:
#        print event.decayName
