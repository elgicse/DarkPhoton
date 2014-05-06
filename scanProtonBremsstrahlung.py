import math
from proton_bremsstrahlung import *
from paraphotonDecay import *
from dark_photon import *
from settings import *
import os

#ngen = 3500

# 0
#mEpsPairs = [[3.e-3, 1.e-3], [9.e-3, 3.e-4], [9.e-3, 2.e-7], [3.e-2, 1.e-4], [8.e-2, 3.e-5], [8.e-2, 7.e-8], [2.e-1, 1.e-5], [2.e-1, 4.e-8], [7.e-1, 1.e-6], [7.e-1, 5.e-8], [9.e-1, 1.e-7]]
# 1
#mEpsPairs = [[3.e-3, 9.e-4], [9.e-3, 2.e-4], [9.e-3, 3.e-7], [3.e-2, 9.e-5], [3.e-2, 1.e-7], [3.e-3, 3.e-7], [8.e-2, 2.e-5], [2.e-1, 9.e-6], [2.e-1, 5.e-8], [7.e-1, 9.e-7], [7.e-1, 6.e-8]]
# 2
#mEpsPairs = [[1.0, 1.e-7], [1.0, 2.e-7], [2.e-1, 1.e-5], [2.e-1, 5.e-8], [8.e-1, 7.e-8], [9.e-1, 4.e-7], [3.e-3, 8.e-4], [3.e-2, 8.e-5], [2.e-1, 8.e-6], [7.e-1, 8.e-7], [1.e-3, 3.e-7]]
# 3
#mEpsPairs = [ [3.e-3, 2.e-7], [1.e-2, 1.e-7], [4.e-2, 8.e-8], [4.e-1, 4.e-7], [6.e-1, 1.e-6], [7.e-2, 2.e-5], [1.e-1,1.5e-5] ]
# 4
#mEpsPairs = [ [2.e-3, 1.e-3], [3.e-3, 7.e-4], [6.e-3, 3.e-4], [3.e-2, 7.e-5], [1.e-1, 1.e-5], [4.e-1, 2.e-6], [2., 2.e-7] ]
# 5
#mEpsPairs = [ [4.e-1, 4.e-8], [1.e-2, 1.5e-7], [1.e-1, 7.e-6], [7.e-1, 1.5e-6], [1.e-2, 1.e-4], [3.e-2, 3.e-5], [3.e-1, 3.e-6] ]
# For the Tracking
mEpsPairs = [[8.e-1, 3.e-7]]


def computeChildrenAcceptance(mass, avgMomentum, pdfTheta, ct, lepton, ngen):
	#print avgMomentum.Pz()
	if not os.path.isfile("out/NTuples/childrenLabFrame_inAcc_%s_m%s.root"%(computeName(lepton), mass)):
		avgBoost = avgMomentum.BoostVector()
		if not os.path.isfile("out/NTuples/childrenRestFrame_%s_m%s.root"%(computeName(lepton), mass)):
			print "Generating rest frame decays..."
			makeNtupleDecayRestFrame(e,mass,ngen)
		f = r.TFile("out/NTuples/childrenRestFrame_%s_m%s.root"%(computeName(lepton), mass), "read")
		dataRF = f.Get("decayRF")
		count2 = 0
		accVol1 = 0
		accVol2 = 0
		res = {}
		dataLF = r.TTree("decayLF","decayLF")
		for (i,children) in enumerate(dataRF):
			if (i%2500)==0:
				print "Boosting RF event %s...."%i
			dataRF.GetEntry(i)
			count2 += 1
			child1 = r.TLorentzVector(dataRF.__getattr__("%s1_Px"%computeName(lepton)),
				dataRF.__getattr__("%s1_Py"%computeName(lepton)),
				dataRF.__getattr__("%s1_Pz"%computeName(lepton)),
				dataRF.__getattr__("%s1_E"%computeName(lepton)))
			child2 = r.TLorentzVector(dataRF.__getattr__("%s2_Px"%computeName(lepton)),
				dataRF.__getattr__("%s2_Py"%computeName(lepton)),
				dataRF.__getattr__("%s2_Pz"%computeName(lepton)),
				dataRF.__getattr__("%s2_E"%computeName(lepton)))
			vtxVol1 = forceDecayInVolume(avgMomentum, pdfTheta, ct, 1)
			vtxVol2 = forceDecayInVolume(avgMomentum, pdfTheta, ct, 2)
			#print vtxVol1.X(), vtxVol1.Y(), vtxVol1.Z()
			#print vtxVol2.X(), vtxVol2.Y(), vtxVol2.Z()
			#avgBoost = avgMomentum.BoostVector()
			#print avgBoost.X(), avgBoost.Y(), avgBoost.Z()
			#print child1.Pz(), ( (+1.)*avgMomentum.Gamma()*avgMomentum.Beta()*child1.E() - child1.Pz()*avgMomentum.Gamma() )
			child1.Boost(avgBoost)
			#print child1.Pz(), child1.Px(), child1.Py(), child1.E(), avgBoost.X(), avgBoost.Y(), avgBoost.Z()
			child2.Boost(avgBoost)
			#print child1.Pz()
			#print avgMomentum.Pz()
			#print inAcceptance(vtxVol1, child1, child2), inAcceptance(vtxVol2, child1, child2)
			a1 = inAcceptance(vtxVol1, child1, child2)
			a2 = inAcceptance(vtxVol2, child1, child2)
			in_acc = 0
			vtx = r.TVector3(0., 0., 0.)
			if a1:
				accVol1 +=1
				in_acc = 1
				vtx = vtxVol1
			if a2:
				accVol2 +=1
				in_acc = 1
				vtx = vtxVol2
			#dataLF.AddVar(avgMomentum.Px(), "avg_A_Px", res)
			#dataLF.AddVar(avgMomentum.Py(), "avg_A_Py", res)
			#dataLF.AddVar(avgMomentum.Pz(), "avg_A_Pz", res)
			#dataLF.AddVar(avgMomentum.E(), "avg_A_E", res)
			#in_acc = int(a1 + a2)
			dataLF.AddVar(vtx.X(), "vtx_x", res)
			dataLF.AddVar(vtx.Y(), "vtx_y", res)
			dataLF.AddVar(vtx.Z(), "vtx_z", res)
			dataLF.AddVar(child1.Px(), "%s1_Px"%computeName(lepton), res)
			dataLF.AddVar(child1.Py(), "%s1_Py"%computeName(lepton), res)
			dataLF.AddVar(child1.Pz(), "%s1_Pz"%computeName(lepton), res)
			dataLF.AddVar(child1.E(),  "%s1_E"%computeName(lepton),  res)
			dataLF.AddVar(child2.Px(), "%s2_Px"%computeName(lepton), res)
			dataLF.AddVar(child2.Py(), "%s2_Py"%computeName(lepton), res)
			dataLF.AddVar(child2.Pz(), "%s2_Pz"%computeName(lepton), res)
			dataLF.AddVar(child2.E(),  "%s2_E"%computeName(lepton),  res)
			dataLF.AddVar(accVol1, "nVol1", res)
			dataLF.AddVar(accVol2, "nVol2", res)
			dataLF.AddVar(in_acc, "in_acceptance", res)
			dataLF.SetDirectory(0)
			dataLF.Fill()
		#print accVol1, accVol2
		accV1 = float(accVol1)/float(count2)
		accV2 = float(accVol2)/float(count2)
		#f.Close()
		ofile = r.TFile("out/NTuples/childrenLabFrame_inAcc_%s_m%s.root"%(computeName(lepton), mass), "recreate")
		dataLF.Write()
		ofile.Close()
	else:
		f = r.TFile("out/NTuples/childrenLabFrame_inAcc_%s_m%s.root"%(computeName(lepton), mass), "read")
		dataLF = f.Get("decayLF")
		count2 = dataLF.GetEntries()
		dataLF.GetEntry(-1)
		accVol1 = dataLF.__getattr__("nVol1")
		accVol2 = dataLF.__getattr__("nVol2")
		accV1 = float(accVol1)/float(count2)
		accV2 = float(accVol2)/float(count2)
		#f.Close()
	return accV1, accV2



def expectedEvents(mEpsPairs, ngen, npp = 1000):
	""" Compute the number of expected events in SHiP for a given array of mass, eps pairs """
	for mass, eps in mEpsPairs:
		# Compute probability that an A' decays in the fiducial volume
		prodFrac = prodRate(mass, eps)
		if not os.path.isfile("out/NTuples/ParaPhoton_eps%s_m%s.root"%(eps, mass)):
			print "Creating four-momenta for eps, mass = ",eps,mass
			create4Momenta(mass, eps, prodFrac, npp)
		f = r.TFile("out/NTuples/ParaPhoton_eps%s_m%s.root"%(eps, mass), "read")
		data = f.Get("newTree")
		count = 0
		probVol1 = 0
		probVol2 = 0
		ct = cTau(mass, eps)
		hx = r.TH1F("hx","hx",50,0.,400.)
		hy = r.TH1F("hy","hy",50,0.,400.)
		hz = r.TH1F("hz","hz",50,0.,400.)
		he = r.TH1F("he","he",50,0.,400.)
		for event in data:
			count += 1
			hx.Fill(event.A_Px)
			hy.Fill(event.A_Py)
			hz.Fill(event.A_Pz)
			he.Fill(event.A_E)
			momentum = r.TLorentzVector(event.A_Px, event.A_Py, event.A_Pz, event.A_E)
			probVol1 += probVtxInVolume(momentum, ct, 1)
			probVol2 += probVtxInVolume(momentum, ct, 2)
		prob1 = probVol1/count
		prob2 = probVol2/count
		#print prob1, prob2
		avgMomentum = r.TLorentzVector(hx.GetMean(), hy.GetMean(), hz.GetMean(), he.GetMean())
		avgBoost = avgMomentum.BoostVector()
		pdfTheta = f.Get("hPDFtheta_eps%s_m%s"%(eps,mass))
		# Compute the probability that the children of an A'
		# decayed in the fiducial volume are detectable
		accV1e, accV2e = computeChildrenAcceptance(mass, avgMomentum, pdfTheta, ct, e, ngen)
		bre = leptonicBranchingRatio(mass, eps, e)
		if mass > 2.*(mmu/1000.):
			accV1mu, accV2mu = computeChildrenAcceptance(mass, avgMomentum, pdfTheta, ct, mu, ngen)
			brmu = leptonicBranchingRatio(mass, eps, mu)
			fracV1 = prob1 * ( bre*accV1e + brmu*accV1mu )
			fracV2 = prob2 * ( bre*accV2e + brmu*accV2mu )
		else:
			fracV1 = prob1 * bre*accV1e
			fracV2 = prob2 * bre*accV2e
		N = protonFlux * prodFrac * (fracV1 + fracV2)
		print "Expected events for m=%s and eps=%s:\t %s"%(mass, eps, N)
		print "protonFlux=%s \t prodFrac=%s \t accV1=%s \t accV2=%s"%(protonFlux, prodFrac, fracV1, fracV2)


