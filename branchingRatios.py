import numpy as np
import ROOT as r
import math
import matplotlib.pyplot as plt
from array import array

from settings import *
from dark_photon import *

masses = np.linspace(0.15,0.75,100).tolist()

e = "e"
mu = "mu"
tau = "tau"

BRe = []
BRmu = []
BRtau = []
BRh = []


def makeBRplot(mDarkPhoton,epsilon):

	hPdfFilename = "out/NTuples/ParaPhoton_eps%s_m%s.root"%(epsilon,mDarkPhoton)
	
	for mass in masses:
		BRe.append(leptonicBranchingRatio(mass, epsilon, e))
		BRmu.append(leptonicBranchingRatio(mass, epsilon, mu))
		BRtau.append(leptonicBranchingRatio(mass, epsilon, tau))
		BRh.append(hadronicBranchingRatio(mass, epsilon))
	
	massesR = array('d',masses)
	BReR = array('d',BRe)
	BRmuR = array('d',BRmu)
	BRtauR = array('d',BRtau)
	BRhR = array('d',BRh)
	
	mgr = r.TMultiGraph("ParaPhotonBranchingRatios_eps%s_m%s"%(epsilon,mDarkPhoton),"ParaPhotonBranchingRatios_eps%s_m%s"%(epsilon,mDarkPhoton))
	gre = r.TGraph(len(massesR),massesR,BReR)
	grmu = r.TGraph(len(massesR),massesR,BRmuR)
	grtau = r.TGraph(len(massesR),massesR,BRtauR)
	grh = r.TGraph(len(massesR),massesR,BRhR)
	gre.SetLineColor(r.kRed)
	gre.SetMarkerColor(r.kRed)
	gre.SetMarkerStyle(20)
	grmu.SetLineColor(r.kBlue)
	grmu.SetMarkerColor(r.kBlue)
	grmu.SetMarkerStyle(20)
	grtau.SetLineColor(r.kGray+2)
	grtau.SetMarkerColor(r.kGray+2)
	grtau.SetMarkerStyle(20)
	grh.SetLineColor(r.kBlack)
	grh.SetMarkerColor(r.kBlack)
	grh.SetMarkerStyle(20)
	mgr.Add(gre)
	mgr.Add(grmu)
	#mgr.Add(grtau)
	mgr.Add(grh)
	
	canvas = r.TCanvas("cParaPhotonBranchingRatios_eps%s_m%s"%(epsilon,mDarkPhoton),"cParaPhotonBranchingRatios_eps%s_m%s"%(epsilon,mDarkPhoton))
	mgr.Draw("alp")
	mgr.GetXaxis().SetTitle("#gamma' mass [GeV]")
	mgr.GetYaxis().SetTitle("Branching ratio")
	mgr.GetYaxis().SetRangeUser(0,1)
	mgr.GetYaxis().SetDecimals()
	
	leg = r.TLegend(0.7, 0.5, 0.9, 0.7)
	leg.SetName("legend")
	leg.SetFillColor(0)
	#leg.SetHeader("#gamma' decay")
	leg.AddEntry(gre, "e^{+}e^{-}", "lp")
	leg.AddEntry(grmu, "#mu^{+}#mu^{-}", "lp")
	leg.AddEntry(grh, "hadrons", "lp")
	leg.Draw()
	
	outfile = r.TFile(hPdfFilename,"update")
	mgr.Write()
	leg.Write()
	canvas.Write()
	outfile.Close()
	
	#plt.ion()
	#plt.scatter(masses,BRe,color='r')
	#plt.scatter(masses,BRmu,color='b')
	#plt.scatter(masses,BRh,color='k')
	#plt.xlabel("A' mass [GeV]")
	#plt.ylabel("Branching ratio")
	
	#plt.show()