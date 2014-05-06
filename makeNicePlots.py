import sys
import os
import ROOT as r

#print sys.argv
mass = sys.argv[2]
eps = sys.argv[1]

r.gROOT.ProcessLine(".x lhcbstyle.C")
filename = "out/NTuples/ParaPhoton_eps%s_m%s.root"%(eps,mass)
infile = r.TFile(filename,"read")
cSaver = []

r.gStyle.SetOptTitle(1)
r.gStyle.SetPadRightMargin(0.08)

c1 = r.TCanvas("chPDF","chPDF")
cSaver.append(c1)
#c1.UseCurrentStyle()
c1.SetTopMargin(0.16)
c1.SetRightMargin(0.16)
pdfhist = infile.Get("hPDF_eps%s_m%s"%(eps,mass))
pdfhist.Draw("COLZ")
title = r.TLatex()
title.SetNDC()
title.SetTextFont(42)
title.SetTextSize(0.06)
title.DrawLatex(0.1,0.9,"PDF for A' production (m_{A'}=%s GeV, #epsilon =%s)"%(mass,eps))
c1.SetLogz()
c1.SetLogx()

c7 = r.TCanvas("chPDFp","chPDFp")
cSaver.append(c7)
pdfhistp = infile.Get("hPDFp_eps%s_m%s"%(eps,mass))
pdfhistp.SetTitle("PDF for A' production (M_{A'}=%s GeV, #epsilon=%s)"%(mass,eps))
pdfhistp.GetXaxis().SetTitle("P_{A'} [GeV]")
pdfhistp.Draw()

c8 = r.TCanvas("chPDFtheta","chPDFtheta")
cSaver.append(c8)
pdfhisttheta = infile.Get("hPDFtheta_eps%s_m%s"%(eps,mass))
pdfhisttheta.SetTitle("PDF for A' production (M_{A'}=%s GeV, #epsilon=%s)"%(mass,eps))
pdfhisttheta.GetXaxis().SetTitle("#theta_{A'} [rad]")
pdfhisttheta.Draw()

r.gStyle.SetOptTitle(0)
r.gStyle.SetOptStat(1110)
tree = infile.Get("newTree")

c2 = r.TCanvas("chvtxz","chvtxz")
cSaver.append(c2)
tree.Draw("A_dec_vtx_z>>hvtxz")
hvtxz = r.gPad.GetPrimitive("hvtxz")
hvtxz.GetXaxis().SetTitle("Paraphoton decay vertex z position [m]")

c3 = r.TCanvas("chvtxx","chvtxx")
cSaver.append(c3)
tree.Draw("A_dec_vtx_x>>hvtxx")
hvtxz = r.gPad.GetPrimitive("hvtxx")
hvtxz.GetXaxis().SetTitle("Paraphoton decay vertex x position [m]")

c4 = r.TCanvas("chtau","chtau")
cSaver.append(c4)
tree.Draw("A_lifetime>>htau")
htau = r.gPad.GetPrimitive("htau")
htau.GetXaxis().SetTitle("Paraphoton lifetime #tau in its rest frame [s]")

c5 = r.TCanvas("chctau","chctau")
cSaver.append(c5)
tree.Draw("A_decLength>>hctau")
hctau = r.gPad.GetPrimitive("hctau")
hctau.GetXaxis().SetTitle("Paraphoton decay length c#tau in its rest frame [m]")

c6 = r.TCanvas("chE","chE")
cSaver.append(c6)
tree.Draw("A_E>>hE")
hE = r.gPad.GetPrimitive("hE")
hE.GetXaxis().SetTitle("Paraphoton energy [GeV]")

c9 = r.TCanvas("chPz","chPz")
cSaver.append(c9)
tree.Draw("A_Pz>>hPz")
hPz = r.gPad.GetPrimitive("hPz")
hPz.GetXaxis().SetTitle("Paraphoton momentum along z (P_{z}) [GeV]")

c10 = r.TCanvas("chPx","chPx")
cSaver.append(c10)
tree.Draw("A_Px>>hPx")
hPx = r.gPad.GetPrimitive("hPx")
hPx.GetXaxis().SetTitle("Paraphoton momentum along x (P_{x}) [GeV]")

c11 = r.TCanvas("cBR","cBR")
cSaver.append(c11)
mgr = infile.Get("ParaPhotonBranchingRatios_eps%s_m%s"%(eps,mass))
leg = infile.Get("legend")
mgr.Draw("alp")
leg.Draw("same")
mgr.GetXaxis().SetLabelFont(42)
mgr.GetXaxis().SetLabelSize(0.05)
mgr.GetXaxis().SetTitleFont(42)
mgr.GetXaxis().SetTitleSize(0.06)
mgr.GetXaxis().SetLabelOffset(0.015)
mgr.GetYaxis().SetLabelFont(42)
mgr.GetYaxis().SetLabelSize(0.05)
mgr.GetYaxis().SetTitleFont(42)
mgr.GetYaxis().SetTitleSize(0.06)
mgr.GetYaxis().SetLabelOffset(0.015)
c11.Modified()
c11.Update()

outputDir = "out/Plots/eps%s_m%s/"%(eps,mass)
if not os.path.exists(outputDir):
    os.makedirs(outputDir)
for canvas in cSaver:
	outprint = outputDir + "eps%s_m%s_"%(eps,mass) + canvas.GetName() + ".pdf"
	canvas.Print(outprint)
	canvas.Close()