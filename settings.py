from __future__ import division
import math

#from functions import *

# dark sector settings
#epsilon = 0.00001
#epsilon = 0.003
#epsilon = 1.0
#mDarkPhoton = 0.07
#mDarkPhoton = 1.0
#mDarkPhoton = 0.05

# useful functions
def energy(p,m):
	""" Compute energy from momentum and mass """
	return math.sqrt(p*p + m*m)

def momentum(E,m):
	""" Compute momentum from energy and mass """
	return math.sqrt(E*E - m*m)

# runtime
res = {}

# proton mass
mProton = 0.938 # GeV/c

# experiment settings
# NuCal1
#protonEnergy = 70. # GeV/c
#protonMomentum = momentum(protonEnergy,mProton)
#protonFlux = 1.71*pow(10.,18.)
#firstVolume = [64., 87., 2.6] # start, end, radius (m)
#secondVolume = [110., 110., 2.6] # start, end, radius (m)
# SHiP
protonEnergy = 400. # GeV/c
protonMomentum = momentum(protonEnergy,mProton)
protonFlux = 2.*pow(10.,20.)
firstVolume = [60., 100., 2.5] # start, end, radius (m)
secondVolume = [110., 150., 2.5] # start, end, radius (m)
v1ThetaMax = firstVolume[2]/firstVolume[0]
v2ThetaMax = secondVolume[2]/secondVolume[0]

# constants
alphaQED = 1./137.
heV = 6.58211928*pow(10.,-16)
hGeV = heV * pow(10.,-9)
c = 3. * pow(10.,8)

# quark masses in MeV
mu = 2.3
md = 4.8
ms = 95.
mc = 1.28 * pow(10.,3.)
mb = 4.18 * pow(10.,3.)
mt = 173. * pow(10.,3.)

# number of colors
ncol = 3.

# quark charges
qu = 2./3.
qd = -1./3.
qs = -1./3.
qc = 2./3.
qb = -1./3.
qt = 2./3.

# lepton masses in MeV
me = 0.511
mmu = 105.7
mtau = 1776.8

# lepton names
e = "e"
mu = "mu"
tau = "tau"

def computeName(lepton):
	if lepton == "e":
		name = "e"
	elif lepton == "mu":
		name = "mu"
	elif lepton == "tau":
		name = "tau"
	else:
		print "ERROR: child lepton undefined. Please select e, mu or tau"
		name = ""
	return name

def computeMass(lepton):
	if lepton == "e":
		ml = me
	elif lepton == "mu":
		ml = mmu
	elif lepton == "tau":
		ml = mtau
	else:
		print "ERROR: child lepton undefined. Please select e, mu or tau"
		ml = -1.
	return ml

# Boundaries for production in meson decays
pi0Start = 0.
pi0Stop = 0.13498
etaStart = 0.13498
etaStop = 0.54785
omegaStart = 0.54785
omegaStop = 0.64767
eta1Start = 0.64767
eta1Stop = 0.95778

# Meson masses
pi0Mass = 0.13498
etaMass = 0.54785
omegaMass = 0.78265
eta1Mass = 0.95778

# Meson decay to gamma branching ratios
def computeMesonBR(mass):
	if pi0Start < mass < pi0Stop:
		# Using pi0 -> gamma gamma
		# From http://pdg.lbl.gov/2013/listings/rpp2013-list-pi-zero.pdf
		br = 98.82
	elif etaStart < mass < etaStop:
		# Using eta -> gamma gamma
		# From http://pdg.lbl.gov/2013/listings/rpp2013-list-eta.pdf
		br = 72.1
	elif omegaStart < mass < omegaStop:
		# Using omega -> pi0 gamma
		# From http://pdg.lbl.gov/2013/listings/rpp2013-list-omega-782.pdf
		br = 8.3
	elif eta1Start < mass < eta1Stop:
		# Using eta' -> gamma gamma
		# From http://pdg.lbl.gov/2013/listings/rpp2013-list-eta-prime-958.pdf
		br = 2.20
	else:
		print "ERROR: invalid A' mass for meson decay production"
		return 0
	return br/100.

def computeScalingFactor(eps, mass):
	phaseSpace = 0
	if pi0Start < mass < pi0Stop:
		motherMass = pi0Mass
	elif etaStart < mass < etaStop:
		motherMass = etaMass
	elif omegaStart < mass < omegaStop:
		motherMass = omegaMass
		phaseSpace = ( 0.5 * pow(motherMass**2. - pi0Mass**2. - mass**2., 2.)
			* math.sqrt( pow((motherMass**2. + pi0Mass**2. - mass**2.), 2.) - 4.*(motherMass**2.)*(pi0Mass**2.) )
			/ pow(motherMass**2. - mass**2., 3.) )
	elif eta1Start < mass < eta1Stop:
		motherMass = eta1Mass
	else:
		print "ERROR: invalid A' mass for meson decay production"
		return 0
	if not phaseSpace:
		phaseSpace = pow((1.- (mass*mass)/(motherMass*motherMass)), 3.)
	br = computeMesonBR(mass)
	factor = 2.*eps*eps*phaseSpace*br
	return factor
