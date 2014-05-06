import numpy as np
import ROOT as r
import math
import matplotlib.pyplot as plt
from array import array

from settings import *
from dark_photon import *

def decayRestFrame(lepton,mDarkPhoton):
	generator = r.TGenPhaseSpace()
	mass = computeMass(lepton)/1000. # in GeV
	if mDarkPhoton < 2.*mass:
		return 0
	else:
		pARestFrame = r.TLorentzVector(0.,0.,0.,mDarkPhoton)
		masses = array('d',[mass,mass])
		generator.SetDecay(pARestFrame, 2, masses)
		generator.Generate()
		pLepton1 = generator.GetDecay(0)
		pLepton2 = generator.GetDecay(1)
		#print pLepton1[0], pLepton1[1], pLepton1[2], pLepton1[3]
		return pLepton1, pLepton2

# Grandchildren to be implemented for a refined search...
def secondDecayRF(lepton):
	generator = r.TGenPhaseSpace()
	if lepton == "mu":
		mass = me
	else:
		print "Error: this is not a muon"
		return 0
	pLrestFrame = r.TLorentzVector(0.,0.,0.,computeMass(lepton))
	masses = array('d',[mass,0.,0.])
	generator.SetDecay(pLrestFrame, 3, masses)
	generator.Generate()
	pGc1 = generator.GetDecay(0)
	pGc2 = generator.GetDecay(1)
	pGc3 = generator.GetDecay(2)

def makeNtupleDecayRestFrame(lepton,mDarkPhoton,nEvents=10000):
	decayRF = r.TTree("decayRF","decayRF")
	#print lepton
	res = {}
	if mDarkPhoton > (2.* computeMass(lepton)/1000.):
		for i in xrange(nEvents):
			if (i%50000)==0:
				print "Making rest frame decay n. %s for A' -> %s%s"%(i,computeName(lepton),computeName(lepton))
			pl1, pl2 = decayRestFrame(lepton,mDarkPhoton)
			#print pl1, pl2
			decayRF.AddVar(pl1.Px(), "%s1_Px"%computeName(lepton), res)
			decayRF.AddVar(pl1.Py(), "%s1_Py"%computeName(lepton), res)
			decayRF.AddVar(pl1.Pz(), "%s1_Pz"%computeName(lepton), res)
			decayRF.AddVar(pl1.E(),  "%s1_E"%computeName(lepton),  res)
			decayRF.AddVar(pl2.Px(), "%s2_Px"%computeName(lepton), res)
			decayRF.AddVar(pl2.Py(), "%s2_Py"%computeName(lepton), res)
			decayRF.AddVar(pl2.Pz(), "%s2_Pz"%computeName(lepton), res)
			decayRF.AddVar(pl2.E(),  "%s2_E"%computeName(lepton),  res)
			decayRF.SetDirectory(0)
			decayRF.Fill()
	outFile = r.TFile("out/NTuples/childrenRestFrame_%s_m%s.root"%(computeName(lepton),mDarkPhoton),"recreate")
	decayRF.Write()
	outFile.Close()


def boostChildren(lepton,mDarkPhoton,epsilon,nMaxParaphotons=0,nMaxChildren=0):
	""" Takes an already existing tree of generated paraphotons, and a tree of decay products in the rest frame.
	For each paraphoton, boosts every decay product with the paraphoton four-momentum. """
	# Input data
	rfFile = r.TFile("out/NTuples/childrenRestFrame_%s_m%s.root"%(computeName(lepton),mDarkPhoton),"read")
	paraphFile = r.TFile("out/NTuples/ParaPhoton_eps%s_m%s.root"%(epsilon,mDarkPhoton),"read")
	rfTree = rfFile.Get("decayRF")
	paraphTree = paraphFile.Get("newTree")
	if not paraphTree:
		paraphTree = paraphFile.Get("newTree;1")
	# Final state tree
	fsTree = r.TTree("finalState","finalState")
	if not nMaxParaphotons:
		nMaxParaphotons = paraphTree.GetEntries()
	if not nMaxChildren:
		nMaxChildren = rfTree.GetEntries()
	nGeneratedEvents = 0
	nEventsInAcceptance = 0
	# Work
	res = {}
	for (j,phot) in enumerate(paraphTree):
		if j > nMaxParaphotons:
			break
		pMother = r.TLorentzVector(phot.A_Px, phot.A_Py, phot.A_Pz, phot.A_E)
		#print "mum: \tTheta = %s \tPz = %s"%(pMother.Theta(),pMother.Pz())
		boostVec = pMother.BoostVector()
		vtx = r.TVector3(phot.A_dec_vtx_x, phot.A_dec_vtx_y, phot.A_dec_vtx_z)
		for (i,childPair) in enumerate(rfTree):
			if i > nMaxChildren:
				break
			rfTree.GetEntry(i)
			nGeneratedEvents += 1
			if nGeneratedEvents%1000000 == 0:
				print "Processing event %s"%nGeneratedEvents
			pChild1 = r.TLorentzVector(rfTree.__getattr__("%s1_Px"%computeName(lepton)),
				rfTree.__getattr__("%s1_Py"%computeName(lepton)),
				rfTree.__getattr__("%s1_Pz"%computeName(lepton)),
				rfTree.__getattr__("%s1_E"%computeName(lepton)),)
			pChild2 = r.TLorentzVector(rfTree.__getattr__("%s2_Px"%computeName(lepton)),
				rfTree.__getattr__("%s2_Py"%computeName(lepton)),
				rfTree.__getattr__("%s2_Pz"%computeName(lepton)),
				rfTree.__getattr__("%s2_E"%computeName(lepton)),)
			#print "%s1: \tTheta = %s \tPz = %s \t\t %s2: \tTheta = %s \tPz = %s"%(computeName(lepton),pChild1.Theta(),pChild1.Pz(),
			#	computeName(lepton),pChild2.Theta(),pChild2.Pz())
			pChild1.Boost(boostVec)
			pChild2.Boost(boostVec)
			#print "%s1b: \tTheta = %s \tPz = %s \t\t %s2b: \tTheta = %s \tPz = %s"%(computeName(lepton),pChild1.Theta(),pChild1.Pz(),
			#	computeName(lepton),pChild2.Theta(),pChild2.Pz())
			inAcc = inAcceptance(vtx, pChild1, pChild2)
			if inAcc:
				nEventsInAcceptance += 1
			# Write final state to tree
			fsTree.AddVar(pMother.Px(), "mother_Px", res)
			fsTree.AddVar(pMother.Py(), "mother_Py", res)
			fsTree.AddVar(pMother.Pz(), "mother_Pz", res)
			fsTree.AddVar(pMother.E(),  "mother_E",  res)
			fsTree.AddVar(vtx.X(), "vtx_x", res)
			fsTree.AddVar(vtx.Y(), "vtx_y", res)
			fsTree.AddVar(vtx.Z(), "vtx_z", res)
			fsTree.AddVar(pChild1.Px(), "child1_Px", res)
			fsTree.AddVar(pChild1.Py(), "child1_Py", res)
			fsTree.AddVar(pChild1.Pz(), "child1_Pz", res)
			fsTree.AddVar(pChild1.E(),  "child1_E",  res)
			fsTree.AddVar(pChild2.Px(), "child2_Px", res)
			fsTree.AddVar(pChild2.Py(), "child2_Py", res)
			fsTree.AddVar(pChild2.Pz(), "child2_Pz", res)
			fsTree.AddVar(pChild2.E(),  "child2_E",  res)
			fsTree.AddVar(inAcc, "in_acceptance", res)
			fsTree.SetDirectory(0)
			fsTree.Fill()
	# Output file
	if nGeneratedEvents:
		fracEventsInAcc = float(nEventsInAcceptance) / float(nGeneratedEvents)
	else:
		fracEventsInAcc = 0
	print "With eps=%s and M=%s:\t\t %s out of %s %s%s events are in SHiP acceptance (%s)"%(epsilon,mDarkPhoton,
		nEventsInAcceptance,nGeneratedEvents,lepton,lepton,fracEventsInAcc)
	labFile = r.TFile("out/NTuples/childrenLabFrame_%s_eps%s_m%s.root"%(computeName(lepton),epsilon,mDarkPhoton),"update")
	fsTree.Write()
	labFile.Close()
	return fracEventsInAcc

def inAcceptance(vtx, pChild1, pChild2):
	# Check if A' -> l+l- vertex is in the first volume
	if (vtx.Z() > firstVolume[0]) and (vtx.Z() < firstVolume[1]):
		if (vtx.X()**2. + vtx.Y()**2.) < firstVolume[2]:
			# Check if child 1 goes through the detector:
			tx1 = pChild1.Px() / pChild1.Pz()
			ty1 = pChild1.Py() / pChild1.Pz()
			endPos1 = r.TVector3()
			endPos1.SetZ(firstVolume[1])
			endPos1.SetX( vtx.X() + tx1*(endPos1.Z() - vtx.Z()) )
			endPos1.SetY( vtx.Y() + ty1*(endPos1.Z() - vtx.Z()) )
			#print "endPos1 ", endPos1.X(), endPos1.Y(), endPos1.Z()
			if (endPos1.X()**2. + endPos1.Y()**2.) < firstVolume[2]:
				# Check if child 2 goes through the detector:
				tx2 = pChild1.Px() / pChild1.Pz()
				ty2 = pChild1.Py() / pChild1.Pz()
				endPos2 = r.TVector3()
				endPos2.SetZ(firstVolume[1])
				endPos2.SetX( vtx.X() + tx2*(endPos2.Z() - vtx.Z()) )
				endPos2.SetY( vtx.Y() + ty2*(endPos2.Z() - vtx.Z()) )
				if (endPos2.X()**2. + endPos2.Y()**2.) < firstVolume[2]:
					return True
	# Check if A' -> l+l- vertex is in the second volume
	elif (vtx.Z() > secondVolume[0]) and (vtx.Z() < secondVolume[1]):
		if (vtx.X()**2. + vtx.Y()**2.) < secondVolume[2]:
			# Check if child 1 goes through the detector:
			tx1 = pChild1.Px() / pChild1.Pz()
			ty1 = pChild1.Py() / pChild1.Pz()
			endPos1 = r.TVector3()
			endPos1.SetZ(secondVolume[1])
			endPos1.SetX( vtx.X() + tx1*(endPos1.Z() - vtx.Z()) )
			endPos1.SetY( vtx.Y() + ty1*(endPos1.Z() - vtx.Z()) )
			if (endPos1.X()**2. + endPos1.Y()**2.) < secondVolume[2]:
				# Check if child 2 goes through the detector:
				tx2 = pChild1.Px() / pChild1.Pz()
				ty2 = pChild1.Py() / pChild1.Pz()
				endPos2 = r.TVector3()
				endPos2.SetZ(secondVolume[1])
				endPos2.SetX( vtx.X() + tx2*(endPos2.Z() - vtx.Z()) )
				endPos2.SetY( vtx.Y() + ty2*(endPos2.Z() - vtx.Z()) )
				if (endPos2.X()**2. + endPos2.Y()**2.) < secondVolume[2]:
					return True
	# Otherwise
	return False