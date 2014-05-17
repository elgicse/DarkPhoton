from __future__ import division
import ROOT as r
import os
import sys
import gc

def makeMothersList(tree):
	mothers = []
	for decay in tree:
		mothers.append(decay.mumId)
	mothers = list(set(mothers))
	return mothers

def makeMassRanges(tree):
	masses = []
	for decay in tree:
		mumMass = decay.mumMass
		childrenMass = 0.
		for i in xrange(len(decay.childrenMass)):
			childrenMass += decay.childrenMass[i]
		masses.append(mumMass - childrenMass)
	masses = list(set(masses))
	masses.sort()
	return masses

def computeRates(tree):
	d = {}
	for (i,decay) in enumerate(tree):
		if not (i%10000):
			print "Reading decay n. %s..."%i
			gc.collect()
		#name = "%s"%(decay.decayName)
		#mumMass = decay.mumMass
		childrenMass = 0.
		for i in xrange(len(decay.childrenMass)):
			childrenMass += decay.childrenMass[i]
		#delta_m = mumMass - childrenMass
		d.setdefault("%s"%(decay.decayName), {"motherMass":decay.mumMass, "massDiff":decay.mumMass - childrenMass, "nPhotons":0})
		d["%s"%(decay.decayName)]["nPhotons"] += 1
		#del name
		#del mumMass
		#del childrenMass
		#del delta_m
	return d

def process(ifile = "out/PythiaData/tree.root"):
	if not os.path.isfile(ifile):
		print "ERROR: please give me a data file :("
		return -1
	f = r.TFile(ifile, "read")
	t = f.Get("events")
	#mothers = makeMothersList(t)
	#masses = makeMassRanges(t)
	d = computeRates(t)
	t.GetEntry(t.GetEntries()-1)
	nGenerated = t.eventNumber
	f.Close()
	return nGenerated, d


# Example sorting
#d1 = sorted(d.items(), key = lambda x: x[1]["nPhoton"], reverse = True)
#d2 = dict(d1)
#masses = []
#for decay in d2.keys():
#	if d2[decay]["massDiff"] > 0:
#    	masses.append(d2[decay]["massDiff"])
