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