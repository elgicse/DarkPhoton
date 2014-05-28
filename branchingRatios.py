import numpy as np
import ROOT as r
import math
import matplotlib.pyplot as plt
from array import array

from settings import *
from dark_photon import *

masses = np.linspace(0.15,0.8,800).tolist()

e = "e"
mu = "mu"
tau = "tau"

BRe = []
BRmu = []
BRtau = []
BRh = []


def makeBRplot(epsilon):

	hPdfFilename = "out/Plots/BRs_eps%s.root"%(epsilon)
	
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
	
	mgr = r.TMultiGraph("ParaPhotonBranchingRatios_eps%s"%(epsilon),"ParaPhotonBranchingRatios_eps%s"%(epsilon))
	gre = r.TGraph(len(massesR),massesR,BReR)
	grmu = r.TGraph(len(massesR),massesR,BRmuR)
	grtau = r.TGraph(len(massesR),massesR,BRtauR)
	grh = r.TGraph(len(massesR),massesR,BRhR)
	gre.SetLineColor(r.kRed)
	gre.SetLineWidth(3)
	gre.SetMarkerColor(r.kRed)
	gre.SetMarkerStyle(7)
	grmu.SetLineColor(r.kBlue)
	grmu.SetLineWidth(3)
	grmu.SetMarkerColor(r.kBlue)
	grmu.SetMarkerStyle(7)
	grtau.SetLineColor(r.kGreen+1)
	grtau.SetLineWidth(3)
	grtau.SetMarkerColor(r.kGreen+1)
	grtau.SetMarkerStyle(7)
	grh.SetLineColor(r.kBlack)
	grh.SetLineWidth(3)
	grh.SetMarkerColor(r.kBlack)
	grh.SetMarkerStyle(7)
	mgr.Add(gre)
	mgr.Add(grmu)
	#mgr.Add(grtau)
	mgr.Add(grh)
	
	canvas = r.TCanvas("cParaPhotonBranchingRatios_eps%s"%(epsilon),"cParaPhotonBranchingRatios_eps%s"%(epsilon))
	r.gStyle.SetTitleFontSize(0.06)
	mgr.SetTitle("#gamma' decay modes")
	mgr.Draw("alp")
	mgr.GetXaxis().SetTitle("#gamma' mass [GeV]")
	mgr.GetYaxis().SetTitle("Branching ratio")
	mgr.GetXaxis().SetTitleSize(0.06)
	mgr.GetXaxis().SetLabelSize(0.05)
	mgr.GetXaxis().SetTitleOffset(0.75)
	mgr.GetYaxis().SetTitleSize(0.06)
	mgr.GetYaxis().SetLabelSize(0.05)
	mgr.GetYaxis().SetTitleOffset(0.82)
	mgr.GetYaxis().SetRangeUser(0,1)
	mgr.GetXaxis().SetRangeUser(0.2,20.)
	mgr.GetYaxis().SetDecimals()
	
	leg = r.TLegend(0.7, 0.6, 0.9, 0.9)
	leg.SetBorderSize(0)
	leg.SetTextSize(0.06)
	leg.SetName("legend")
	leg.SetFillColor(0)
	#leg.SetHeader("#gamma' decay")
	leg.AddEntry(gre, "e^{+}e^{-}", "lp")
	leg.AddEntry(grmu, "#mu^{+}#mu^{-}", "lp")
	#leg.AddEntry(grtau, "#tau^{+}#tau^{-}", "lp")
	leg.AddEntry(grh, "hadrons", "lp")
	leg.Draw()
	
	outfile = r.TFile(hPdfFilename,"recreate")
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