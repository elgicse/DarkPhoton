import ROOT as r
from array import array


if __name__ == '__main__':

	current_mass = []
	current_eps = []
	
	with open("current_limits_all.csv","r") as f:
		for line in f:
			line = line.split(",")
			current_mass.append(float(line[0]))
			current_eps.append(float(line[1]))
	
	current_r_mass = array('f', current_mass)
	current_r_eps = array('f', current_eps)
	
	gr = r.TGraph(len(current_r_mass), current_r_mass, current_r_eps)
	c3 = r.TCanvas()
	c3.cd()
	gr.SetLineWidth(3504)
	gr.SetFillStyle(3002)
	gr.SetLineColor(r.kGray+2)
	gr.SetTitle("Current limits on hidden photons")
	gr.Draw("alp")
	gr.GetXaxis().SetTitle(r"m_{#gamma'} (GeV)")
	gr.GetYaxis().SetTitle(r"#chi")
	gr.GetXaxis().SetTitleSize(0.05)
	gr.GetYaxis().SetTitleSize(0.05)
	gr.GetXaxis().SetTitleOffset(0.90)
	gr.GetYaxis().SetTitleOffset(0.90)
	c3.SetGrid()
	c3.SetLogx()
	c3.SetLogy()