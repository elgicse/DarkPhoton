import numpy as np
import ROOT as r
import math
from scipy.integrate import quad, dblquad

from settings import *
from dark_photon import *

def zeta(p,theta):
	""" Fraction of the proton momentum carried away by the paraphoton in the beam direction """
	#return 1. / math.sqrt(theta*theta + 1.)
	return p / (protonMomentum * math.sqrt(theta*theta + 1.))

def pTransverse(p,theta):
	""" Paraphoton transverse momentum in the lab frame """
	#return p*theta*zeta(p,theta)
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
	p1 = a2*pow(math.log(s/a3),2.) # natural logarithm, right?
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
	#pPar = p*zeta(p,theta)
	#return 2.*math.fabs(theta)*pPar*pPar
	z2 = pow(zeta(p,theta),2.)
	return 2.*theta*z2*protonMomentum*protonMomentum

def dZdP(p,theta):
	""" Jacobian z->p """
	#return p / (zeta(p,theta)*protonMomentum*protonMomentum)
	#return 1. / (zeta(p,theta)*p)
	return 1./( protonMomentum* math.sqrt(theta*theta+1.) )

def dNdPdTheta(p,theta,mDarkPhoton,epsilon):
	""" Differential A' rate per p.o.t. as a function of P and theta """
	diffRate = dNdZdPtSquare(p,mDarkPhoton,theta,epsilon) * dPt2dTheta(p,theta) * dZdP(p,theta)
	return math.fabs(diffRate) # integrating in (-pi, pi)...

def prodRate(mDarkPhoton,epsilon):
	""" dNdPdTheta integrated over p and theta """
	integral = dblquad( dNdPdTheta, # integrand
						-1*math.pi, 1.*math.pi, # theta boundaries (2nd argument of integrand)
						lambda x: 0., lambda x: protonMomentum, # p boundaries (1st argument of integrand)
						args=(mDarkPhoton, epsilon) ) # extra parameters to pass to integrand
	return integral[0]

# total production rate of A'
norm = prodRate(mDarkPhoton,epsilon)
# number of A' produced
numDarkPhotons = int(math.floor(norm*protonFlux))

print
print "Epsilon \t %s"%epsilon
print "mDarkPhoton \t %s"%mDarkPhoton
print "A' production rate per p.o.t: \t %.8g"%norm
print "Number of A' produced in SHiP: \t %.8g"%numDarkPhotons

def normalisedProductionPDF(p,theta,mDarkPhoton,epsilon,norm):
	""" Probability density function for A' production in SHIP """
	return (1./norm)*dNdPdTheta(p,theta,mDarkPhoton,epsilon)

def hProdPDF(mDarkPhoton,epsilon,norm):
	""" Histogram of the PDF for A' production in SHIP """
	angles = np.linspace(-1*math.pi,1.*math.pi,360,endpoint=False).tolist()
	anglestep = math.pi/180.
	momentumStep = 0.05 # GeV
	ndiv = int(math.floor(protonMomentum/momentumStep))
	momenta = np.linspace(momentumStep,protonMomentum,ndiv,endpoint=False).tolist()
	hPDF = r.TH2F("hPDF_eps%s_m%s"%(epsilon,mDarkPhoton),"hPDF_eps%s_m%s"%(epsilon,mDarkPhoton),
		ndiv,0.5*momentumStep,protonMomentum-0.5*momentumStep,
		360,-1*math.pi-0.5*anglestep,1.*math.pi-0.5*anglestep)
	hPDF.SetTitle("PDF for A' production (m_{A'}=%s GeV, #epsilon =%s)"%(mDarkPhoton,epsilon))
	hPDF.GetXaxis().SetTitle("P_{A'} [GeV]")
	hPDF.GetYaxis().SetTitle("#theta_{A'} [rad]")
	hPDFtheta = r.TH1F("hPDFtheta_eps%s_m%s"%(epsilon,mDarkPhoton),
		"hPDFtheta_eps%s_m%s"%(epsilon,mDarkPhoton),
		360,-1*math.pi-0.5*anglestep,1.*math.pi-0.5*anglestep)
	hPDFp = r.TH1F("hPDFp_eps%s_m%s"%(epsilon,mDarkPhoton),
		"hPDFp_eps%s_m%s"%(epsilon,mDarkPhoton),
		ndiv,0.5*momentumStep,protonMomentum-0.5*momentumStep)
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
	hPDF.Write()
	hPDFp.Write()
	hPDFtheta.Write()
	outfile.Close()
	return hPDF

def create4Momenta(mDarkPhoton,epsilon,norm,nEvents=10000):
	""" Create TTree containing the production PDF and the four-momenta and decay vertices of 1000 A' produced in the SHIP beam dump """
	print "Producing PDF, be patient..."
	pdf = hProdPDF(mDarkPhoton,epsilon,norm)
	#proton4mom = r.TLorentzVector(0., 0., protonMomentum, protonEnergy)
	#boostVect = proton4mom.BoostVector()
	origin = r.TVector3(0.,0.,0.)
	newTree = r.TTree("newTree","newTree")
	res = {}
	for i in xrange(nEvents):
		if i%1000 == 0:
			print "Processing paraphoton %s..."%i
		p, theta = r.Double(), r.Double()
		pdf.GetRandom2(p, theta)
		#print "Momentum %s \t\t Angle %s"%(p,theta)
		vec = r.TVector3()
		vec.SetMagThetaPhi(p, theta, 0) # Phi is arbitrary
		pParaphoton = r.TLorentzVector()
		pParaphoton.SetE(energy(p, mDarkPhoton))
		pParaphoton.SetVect(vec)
		#pParaphoton.Boost(boostVect)

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
