import numpy as np
import ROOT as r
import math
import os
from scipy.integrate import quad, dblquad

from settings import *
from dark_photon import *
from paraphotonDecay import *

def zeta(p,theta):
	""" Fraction of the proton momentum carried away by the paraphoton in the beam direction """
	return p / (protonMomentum * math.sqrt(theta*theta + 1.))

def pTransverse(p,theta):
	""" Paraphoton transverse momentum in the lab frame """
	return protonMomentum*theta*zeta(p,theta)

def ptSquare(p,theta):
	""" Square paraphoton transverse momentum in the lab frame """
	return pow(pTransverse(p,theta), 2.)

def H(p,theta,mDarkPhoton):
	""" A kinematic term """
	return ptSquare(p,theta) + (1.-zeta(p,theta))*mDarkPhoton*mDarkPhoton + pow(zeta(p,theta),2.)*mProton*mProton

def wba(p,theta,mDarkPhoton,epsilon):
	""" Cross section weighting function in the Fermi-Weizsaeker-Williams approximation """
	const = epsilon*epsilon*alphaQED / (2.*math.pi*H(p,theta,mDarkPhoton))

	h2 = pow(H(p,theta,mDarkPhoton),2.)
	oneMinusZSquare = pow(1.-zeta(p,theta),2.)
	mp2 = mProton*mProton
	mA2 = mDarkPhoton*mDarkPhoton

	p1 = (1. + oneMinusZSquare) / zeta(p,theta)
	p2 = ( 2. * zeta(p,theta) * (1.-zeta(p,theta)) * ( (2.*mp2 + mA2)/ H(p,theta,mDarkPhoton) 
			- pow(zeta(p,theta),2.)*2.*mp2*mp2/h2 ) )
	p3 = 2.*zeta(p,theta)*(1.-zeta(p,theta))*(zeta(p,theta)+oneMinusZSquare)*mp2*mA2/h2
	p4 = 2.*zeta(p,theta)*oneMinusZSquare*mA2*mA2/h2
	return const*(p1+p2+p3+p4)

def sigma(s): # s in GeV^2 ---> sigma in mb
	""" Parametrisation of sigma(s) """
	a1 = 35.45
	a2 = 0.308
	a3 = 28.94
	a4 = 33.34
	a5 = 0.545
	a6 = 0.458
	p1 = a2*pow(math.log(s/a3),2.) 
	p2 = a4*pow((1./s),a5)
	p3 = a4*pow((1./s),a6)
	return a1 + p1 - p2 + p3

def es(p,mDarkPhoton):
	""" s(p,mA) """
	return 2.*mProton*(energy(protonMomentum,mProton)-energy(p,mDarkPhoton))

def sigmaRatio(p,mDarkPhoton):
	""" sigma(s') / sigma(s) """
	return sigma(es(p,mDarkPhoton)) / sigma(2.*mProton*energy(protonMomentum,mProton))

def dNdZdPtSquare(p,mDarkPhoton,theta,epsilon):
	""" Differential A' rate per p.o.t. as a function of Z and Pt^2 """
	return sigmaRatio(p,mDarkPhoton)*wba(p,theta,mDarkPhoton,epsilon)

def dPt2dTheta(p,theta):
	""" Jacobian Pt^2->theta """
	z2 = pow(zeta(p,theta),2.)
	return 2.*theta*z2*protonMomentum*protonMomentum

def dZdP(p,theta):
	""" Jacobian z->p """
	return 1./( protonMomentum* math.sqrt(theta*theta+1.) )

def dNdPdTheta(p,theta,mDarkPhoton,epsilon):
	""" Differential A' rate per p.o.t. as a function of P and theta """
	diffRate = dNdZdPtSquare(p,mDarkPhoton,theta,epsilon) * dPt2dTheta(p,theta) * dZdP(p,theta)
	return math.fabs(diffRate) # integrating in (-pi, pi)...

def pMin(mDarkPhoton):
	return max(0.14*protonMomentum, mDarkPhoton)

def pMax(mDarkPhoton):
	#return min(0.86*protonMomentum, math.sqrt( (energy(protonMomentum,mProton)**2. - mDarkPhoton**2.) - mDarkPhoton**2.))
	return math.sqrt( (energy(protonMomentum,mProton)**2. - mDarkPhoton**2.) - mDarkPhoton**2.)

def prodRate(mDarkPhoton,epsilon, tmin = -0.5*math.pi, tmax = 0.5*math.pi):
	""" dNdPdTheta integrated over p and theta """
	integral = dblquad( dNdPdTheta, # integrand
						tmin, tmax, # theta boundaries (2nd argument of integrand)
						lambda x: pMin(mDarkPhoton), lambda x: pMax(mDarkPhoton), # p boundaries (1st argument of integrand)
						args=(mDarkPhoton, epsilon) ) # extra parameters to pass to integrand
	return integral[0]

# total production rate of A'
#norm = prodRate(mDarkPhoton,epsilon)
# number of A' produced
#numDarkPhotons = int(math.floor(norm*protonFlux))
#
#print
#print "Epsilon \t %s"%epsilon
#print "mDarkPhoton \t %s"%mDarkPhoton
#print "A' production rate per p.o.t: \t %.8g"%norm
#print "Number of A' produced in SHiP: \t %.8g"%numDarkPhotons

def normalisedProductionPDF(p,theta,mDarkPhoton,epsilon,norm):
	""" Probability density function for A' production in SHIP """
	return (1./norm)*dNdPdTheta(p,theta,mDarkPhoton,epsilon)

def hProdPDF(mDarkPhoton,epsilon,norm,binsp,binstheta,tmin = -0.5*math.pi, tmax = 0.5*math.pi):
	""" Histogram of the PDF for A' production in SHIP """
	angles = np.linspace(tmin,tmax,binstheta).tolist()
	anglestep = 2.*(tmax - tmin)/binstheta
	momentumStep = (pMax(mDarkPhoton)-pMin(mDarkPhoton))/(binsp-1)
	momenta = np.linspace(pMin(mDarkPhoton),pMax(mDarkPhoton),binsp,endpoint=False).tolist()
	hPDF = r.TH2F("hPDF_eps%s_m%s"%(epsilon,mDarkPhoton) ,"hPDF_eps%s_m%s"%(epsilon,mDarkPhoton),
		binsp,pMin(mDarkPhoton)-0.5*momentumStep,pMax(mDarkPhoton)-0.5*momentumStep,
		binstheta,tmin-0.5*anglestep,tmax-0.5*anglestep)
	hPDF.SetTitle("PDF for A' production (m_{A'}=%s GeV, #epsilon =%s)"%(mDarkPhoton,epsilon))
	hPDF.GetXaxis().SetTitle("P_{A'} [GeV]")
	hPDF.GetYaxis().SetTitle("#theta_{A'} [rad]")
	hPDFtheta = r.TH1F("hPDFtheta_eps%s_m%s"%(epsilon,mDarkPhoton),
		"hPDFtheta_eps%s_m%s"%(epsilon,mDarkPhoton),
		binstheta,tmin-0.5*anglestep,tmax-0.5*anglestep)
	hPDFp = r.TH1F("hPDFp_eps%s_m%s"%(epsilon,mDarkPhoton),
		"hPDFp_eps%s_m%s"%(epsilon,mDarkPhoton),
		binsp,pMin(mDarkPhoton)-0.5*momentumStep,pMax(mDarkPhoton)-0.5*momentumStep)
	hPDFp.GetXaxis().SetTitle("P_{A'} [GeV]")
	hPDFtheta.GetXaxis().SetTitle("#theta_{A'} [rad]")
	for theta in angles:
		for p in momenta:
			w = normalisedProductionPDF(p,theta,mDarkPhoton,epsilon,norm)
			hPDF.Fill(p,theta,w)
			hPDFtheta.Fill(theta,w)
			hPDFp.Fill(p,w)
	hPdfFilename = "out/NTuples/ParaPhoton_eps%s_m%s.root"%(epsilon,mDarkPhoton)
	outfile = r.TFile(hPdfFilename,"recreate")
	weight = hPDF.Integral("width")
	hPDF.Scale(1./weight)
	hPDF.Write()
	hPDFp.Write()
	hPDFtheta.Write()
	outfile.Close()
	del angles
	del momenta
	return hPDF

def create4Momenta(mDarkPhoton,epsilon,norm,nEvents=10000,tmin = -0.5*math.pi, tmax = 0.5*math.pi):
	""" Create TTree containing the production PDF and the four-momenta and decay vertices of 1000 A' produced in the SHIP beam dump """
	print "Producing PDF, be patient..."
	pdf = hProdPDF(mDarkPhoton,epsilon,norm,binsp, binstheta,tmin, tmax)
	origin = r.TVector3(0.,0.,0.)
	newTree = r.TTree("newTree","newTree")
	res = {}
	for i in xrange(nEvents):
		if i%1000 == 0:
			print "Processing paraphoton %s..."%i
		p, theta = r.Double(), r.Double()
		pdf.GetRandom2(p, theta)
		vec = r.TVector3()
		vec.SetMagThetaPhi(p, theta, 0) # Phi is arbitrary
		pParaphoton = r.TLorentzVector()
		pParaphoton.SetE(energy(p, mDarkPhoton))
		pParaphoton.SetVect(vec)

		decayLength = cTau(mDarkPhoton, epsilon)
		tau = lifetime(mDarkPhoton, epsilon)
		vtx = GetDecayPoint(origin, pParaphoton, tau)

		newTree.AddVar(pParaphoton.Px(), "A_Px", res)
		newTree.AddVar(pParaphoton.Py(), "A_Py", res)
		newTree.AddVar(pParaphoton.Pz(), "A_Pz", res)
		newTree.AddVar(pParaphoton.E(), "A_E", res)
		newTree.AddVar(decayLength, "A_decLength", res)
		newTree.AddVar(tau, "A_lifetime", res)
		newTree.AddVar(vtx.X(), "A_dec_vtx_x", res)
		newTree.AddVar(vtx.Y(), "A_dec_vtx_y", res)
		newTree.AddVar(vtx.Z(), "A_dec_vtx_z", res)
		newTree.AddVar(norm, "production_rate", res)

		newTree.SetDirectory(0)
		newTree.Fill()

	hPdfFilename = "out/NTuples/ParaPhoton_eps%s_m%s.root"%(epsilon,mDarkPhoton)
	outfile = r.TFile(hPdfFilename,"update")
	newTree.Write()
	outfile.Close()


def scanPDF(mass, eps, mesonDecay=False):
	outData = []
	ct = cTau(mass, eps)
	# Use the right file
	if not mesonDecay:
		f = r.TFile("out/NTuples/ParaPhoton_eps%s_m%s.root"%(eps,mass),
			"update")
		hpdf = f.Get("hPDF_eps%s_m%s"%(eps,mass))
	elif mesonDecay:
		f = r.TFile("out/PythiaData/mesonDecays_%s.root"%mass,
			"update")
		hpdforig = f.Get("hPDF_m%s"%mass)
		try:
			hpdf = hpdforig.Clone("hpdf")
		except ReferenceError:
			print mass, eps
			hpdforig = f.Get("hPDF_m%s;1"%mass)
			hpdf = hpdforig.Clone("hpdf")
		integr = hpdf.Integral("width")
		hpdf.Scale(1./integr)
	# Get the binning of the original PDF
	DeltaTheta = (hpdf.GetYaxis().GetBinCenter(1)
		-hpdf.GetYaxis().GetBinCenter(0))
	DeltaP = (hpdf.GetXaxis().GetBinCenter(1)
		-hpdf.GetXaxis().GetBinCenter(0))
	ThetaMax = hpdf.GetYaxis().GetXmax()
	ThetaMin = hpdf.GetYaxis().GetXmin()
	RangeTheta = ThetaMax - ThetaMin
	P_Max = hpdf.GetXaxis().GetXmax()
	P_Min = hpdf.GetXaxis().GetXmin()
	RangeP = P_Max - P_Min 
	binRelSize = (DeltaP*DeltaTheta)#/(RangeTheta/RangeP)
	nMom = hpdf.GetXaxis().GetNbins()
	nTheta = hpdf.GetYaxis().GetNbins()
	fourMom = r.TLorentzVector()
	vec = r.TVector3()
	# Make PDFs rescaled to the acceptance of the experiment
	hPdfAcc1 = r.TH2F("hPDFinAcc1_eps%s_m%s"%(eps,mass),
		"hPDFinAcc1_eps%s_m%s"%(eps,mass),
		nMom,P_Min,P_Max,nTheta,ThetaMin,ThetaMax)
	hPdfAcc2 = r.TH2F("hPDFinAcc2_eps%s_m%s"%(eps,mass),
		"hPDFinAcc2_eps%s_m%s"%(eps,mass),
		nMom,P_Min,P_Max,nTheta,ThetaMin,ThetaMax)
	index = 0
	valAcc1 = 0.
	valAcc2 = 0.
	wtot = 0.
	naccgeo = 0
	for p in xrange(nMom):
		for th in xrange(nTheta):
			index += 1
			#if (index%1000 == 0):
				#print "bin %s..."%index
			weight = hpdf.GetBinContent(p, th)
			if weight > 0.:
				mom = hpdf.GetXaxis().GetBinCenter(p)
				angle = hpdf.GetYaxis().GetBinCenter(th)
				binWeight = weight*binRelSize
				wtot = wtot + binWeight
				vec.SetMagThetaPhi(mom, angle, 0.)
				fourMom.SetE(energy(mom, mass))
				fourMom.SetVect(vec)
				gamma = fourMom.Gamma()
				px = fourMom.Px()
				pz = fourMom.Pz()
				accGeo1 = GeometricAcceptance(px, pz, 1)
				accGeo2 = GeometricAcceptance(px, pz, 2)
				accLifetimeVol1 = probVtxInVolume(fourMom, ct, 1, gamma)
				accLifetimeVol2 = probVtxInVolume(fourMom, ct, 2, gamma)
				acc1 = binWeight * accGeo1 * accLifetimeVol1
				valAcc1 = valAcc1 + acc1
				acc2 = binWeight * accGeo2 * accLifetimeVol2
				valAcc2 = valAcc2 + acc2
				outData.append([mom, angle, gamma, ct, binWeight, accGeo1,
					accGeo2, accLifetimeVol1, accLifetimeVol2])
				hPdfAcc1.Fill(mom, angle, acc1)
				hPdfAcc2.Fill(mom, angle, acc2)
	if valAcc1 > 1.e-20:
		totWeight1 = hPdfAcc1.Integral("width")
		normalization1 = 1./totWeight1
		hPdfAcc1.Scale(normalization1)
		hPdfAcc1.Write("",5)
	else:
		valAcc1 = 0.
	if valAcc2 > 1.e-20:
		totWeight2 = hPdfAcc2.Integral("width")
		normalization2 = 1./totWeight2
		hPdfAcc2.Scale(normalization2)
		hPdfAcc2.Write("",5)
	else:
		valAcc2 = 0.
	f.Close()
	return valAcc1, valAcc2, outData


def makeAcceptancePdf(mass, eps, binsp, binstheta, mesonDecay):
	tmax = v1ThetaMax
	tmin = 0.#(-1.)*tmax
	if not mesonDecay:
		norm = prodRate(mass, eps, tmin, tmax)
	else:
		norm = prodRateFromMesons(mass)
	if (not mesonDecay) and (not os.path.isfile("out/NTuples/ParaPhoton_eps%s_m%s.root"%(eps,mass))):
		hProdPDF(mass, eps, norm, binsp, binstheta, tmin, tmax)
	valAcc1, valAcc2, outData = scanPDF(mass, eps, mesonDecay)
	return norm, valAcc1, valAcc2, outData


def computeNEvents(mass, eps, mesonDecay=False, binsp=90, binstheta=80):
	outData = makeAcceptancePdf(mass, eps, binsp, binstheta, mesonDecay)
	prodFrac = outData[0]
	# Meson decay ntuples also produce A's out of the acceptances
	if mesonDecay:
		# This is what we do with p -> p+A',
		# we integrate the PDF in the acceptance of volume 1
		prodFrac = prodFrac * geomAcceptance(mass, 1)
	# P(vtx in volume)
	prob1 = outData[1]
	prob2 = outData[2]
	if prob1 or prob2:
		makeNtupleDecayRestFrame(e,mass,200)
	if prob1:
		acc1e = boostChildrenInAcceptance(e,mass,eps,1,mesonDecay,200)
	else:
		acc1e = 0.
	if prob2:
		acc2e = boostChildrenInAcceptance(e,mass,eps,2,mesonDecay,200)
	else:
		acc2e = 0.
	bre = leptonicBranchingRatio(mass, eps, e)
	if mass > 2.*mmu/1000.:
		if prob1 or prob2:
			makeNtupleDecayRestFrame(mu,mass,200)
		if prob1:
			acc1mu = boostChildrenInAcceptance(mu,mass,eps,1,mesonDecay,200)
		else:
			acc1mu = 0.
		if prob2:
			acc2mu = boostChildrenInAcceptance(mu,mass,eps,2,mesonDecay,200)
		else:
			acc2mu = 0.
		brmu = leptonicBranchingRatio(mass, eps, mu)
		fracV1 = prob1 * ( bre*acc1e + brmu*acc1mu )
		fracV2 = prob2 * ( bre*acc2e + brmu*acc2mu )
	else:
		fracV1 = prob1 * bre*acc1e
		fracV2 = prob2 * bre*acc2e
		brmu, acc1mu, acc2mu = 0., 0., 0.
	expectedEvents = protonFlux * prodFrac * (fracV1 + fracV2)
	if mesonDecay:
		factor = computeScalingFactor(eps, mass)
		expectedEvents = expectedEvents * factor
		outFilePath = "out/TextData/sensitivityScan-MesonDecays.txt"
	else:
		outFilePath = "out/TextData/sensitivityScan-FWapprox.txt"
	with open(outFilePath,"a") as ofile:
	#with open("out/TextData/sensitivityScanNuCal1.txt","a") as ofile:
		try:
			if mesonDecay:
				ofile.write("%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s\n"%(mass, eps, prodFrac, prob1, prob2, bre, brmu, acc1e, acc2e, acc1mu, acc2mu, factor, expectedEvents ))
			else:
				ofile.write("%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s\n"%(mass, eps, prodFrac, prob1, prob2, bre, brmu, acc1e, acc2e, acc1mu, acc2mu, expectedEvents ))
		except KeyboardInterrupt:
			pass
	return expectedEvents









def printAcceptancePdfFile(mass, eps, binsp, binstheta):
	outData = makeAcceptancePdf(mass, eps, binsp, binstheta)
	with open("out/TextData/ParaphotonInAcceptancePDF.dat","w") as ofile:
		ofile.write("P \t Theta \t Gamma \t cTau \t binWeight \t accGeo1 \t accGeo2 \t accLifetimeVol1 \t accLifetimeVol2\n")
		for datalist in outData[3]:
			ofile.write("%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t "%(datalist[0],datalist[1],datalist[2],datalist[3],datalist[4],datalist[5],datalist[6],datalist[7],datalist[8]))


