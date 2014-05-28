import math

from settings import *
from functions import *

def leptonicDecayWidth(mDarkPhoton, epsilon, lepton): # mDarkPhoton in GeV
	""" Dark photon decay width into leptons, in GeV (input short name of lepton family) """
	ml = computeMass(lepton)
	ml = ml/1000. # all in GeV
	constant = (1./3.) * alphaQED * mDarkPhoton * pow(epsilon, 2.)
	if 2.*ml < mDarkPhoton:
		rad = math.sqrt( 1. - (4.*ml*ml)/(mDarkPhoton*mDarkPhoton) )
	else:
		rad = 0.
	par = 1. + (2.*ml*ml)/(mDarkPhoton*mDarkPhoton)
	return math.fabs(constant*rad*par)

def leptonicBranchingRatio(mDarkPhoton, epsilon, lepton):
	return leptonicDecayWidth(mDarkPhoton, epsilon, lepton) / totalDecayWidth(mDarkPhoton, epsilon)
	#if lepton == e:
	#	otherlepton = mu
	#	br = 1. - ((leptonicDecayWidth(mDarkPhoton, epsilon, otherlepton) + hadronicDecayWidth(mDarkPhoton, epsilon))/totalDecayWidth(mDarkPhoton, epsilon))
	#else: 
	#	br = leptonicDecayWidth(mDarkPhoton, epsilon, lepton) / totalDecayWidth(mDarkPhoton, epsilon)
	#return br

def hadronicDecayWidth(mDarkPhoton, epsilon):
	""" Dark photon decay into hadrons """
	constant = (1./3.) * alphaQED * mDarkPhoton * pow(epsilon, 2.)
	#if Ree_const in vars():
	#	return constant*Ree_const
	#else:
	#MeVmass = mDarkPhoton * 1000.
	return constant*Ree_interp(mDarkPhoton)

def hadronicBranchingRatio(mDarkPhoton, epsilon):
	return hadronicDecayWidth(mDarkPhoton, epsilon) / totalDecayWidth(mDarkPhoton, epsilon)

def totalDecayWidth(mDarkPhoton, epsilon): # mDarkPhoton in GeV
	""" Total decay width in GeV """
	#return hGeV*c / cTau(mDarkPhoton, epsilon)
	tdw = (leptonicDecayWidth(mDarkPhoton, epsilon, e)
		+ leptonicDecayWidth(mDarkPhoton, epsilon, mu)
		+ leptonicDecayWidth(mDarkPhoton, epsilon, tau)
		+ hadronicDecayWidth(mDarkPhoton, epsilon))
	return tdw
	#if Ree_const in vars():
	#	return 3./((1.+Ree_const)*alphaQED*mDarkPhoton*pow(epsilon,2.))
	#else:
	#	return 3./((1.+Ree_interp(mDarkPhoton))*alphaQED*mDarkPhoton*pow(epsilon,2.))

def cTau(mDarkPhoton, epsilon): # decay length in meters, dark photon mass in GeV
	""" Dark Photon lifetime """
	#MeVmass = mDarkPhoton * 1000.
	#if Ree_const in vars():
	#	p1 = 0.8 / (1. + Ree_const)
	#else:
	p1 = 0.8 / (1. + Ree_interp(mDarkPhoton))
	p2 = pow((pow(10., -6.) / epsilon),2.)
	p3 = 100. / (mDarkPhoton*1000.)
	#print "ctau ",p1*p2*p3, 
	return p1*p2*p3 #GeV/MeV conversion

def lifetime(mDarkPhoton, epsilon):
	c = 3.*pow(10.,8.)
	return cTau(mDarkPhoton, epsilon)/c

def lifetime100ms(mass):
	rad = 0.8/(3.*mass*(1+Ree_interp(mass)))
	return math.sqrt(rad)*pow(10.,-10.5)


#def Ree(s): # s in MeV
#	""" sigma(e+e- -> hadrons) / sigma(e+e- -> mu+mu-) """
#	if s > 2.*mt:
#		ratio = ncol*(qu*qu + qd*qd + qc*qc + qs*qs + qb*qb + qt*qt)
#	elif s > 2.*mb:
#		ratio = ncol*(qu*qu + qd*qd + qc*qc + qs*qs + qb*qb)
#	elif s > 2.*mc:
#		ratio = ncol*(qu*qu + qd*qd + qc*qc + qs*qs)
#	elif s > 2.*ms:
#		ratio = ncol*(qu*qu + qd*qd + qs*qs)
#	elif s > 2.*md:
#		ratio = ncol*(qu*qu + qd*qd)
#	elif s > 2.*mu:
#		ratio = ncol*(qu*qu)
#	else:
#		ratio = 0

def Ree_interp(s): # s in GeV
	""" Using PDG values for sigma(e+e- -> hadrons) / sigma(e+e- -> mu+mu-) """
	# Da http://pdg.lbl.gov/2012/hadronic-xsections/hadron.html#miscplots
	#ecm = math.sqrt(s)
	ecm = s
	if ecm>=dataEcm[0]:
		result = float(PdgR(ecm))
	else:
		result=0
	return result
