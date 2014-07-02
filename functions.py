from __future__ import division
import math
import os
import ROOT as r
import numpy as np
from array import array
from scipy import interpolate 
import matplotlib.pyplot as plt

from settings import *


def __AddVar__(self, var, name, res):
    has_branche = res.has_key(name)
            
    if has_branche:
        res[name][0] = float(var)
    else:
        comp = array('f', [0])
        self.Branch(name, comp, name+'/F')
        comp[0]= float(var)
        aux = {name: comp}
        res.update(aux)
    return res

r.TTree.AddVar = __AddVar__

def GetDecayPoint(Origin, Momentum, lifetime):
    c = 3.*pow(10.,8.)
    v = Momentum.Beta()*c  
    DT =r.gRandom.Exp(lifetime)*Momentum.Gamma()
    DL = DT*v
    Direction = Momentum.Vect().Unit()
    EndVertex = r.TVector3(Momentum.Vect().Unit()[0]*DL, Momentum.Vect().Unit()[1]*DL, Momentum.Vect().Unit()[2]*DL) 
    EndVertex = Origin+EndVertex
    return EndVertex


def probVtxInVolume(momentum, ct, volume, gamma=0):
	if volume == 1:
		vol = firstVolume
	elif volume == 2:
		vol = secondVolume
	else:
		print "ERROR: select decay volume 1 or 2"
		return 0
	if not gamma:
		gamma = momentum.Gamma()
	costheta = np.fabs(momentum.CosTheta())
	tantheta = np.tan(momentum.Theta())
	start = vol[0]
	end = vol[1]
	rad = vol[2]
	stop = min(rad/tantheta, end)
	if stop < start:
		return 0.
	esp1 = (-1.) * (start/costheta) / (gamma*ct)
	esp2 = (-1.) * (stop/costheta) / (gamma*ct)
	#np.seterr(all='raise')
	try:
		result = np.nan_to_num(np.fabs( np.exp(esp1) - np.exp(esp2) ))
	except (ValueError, FloatingPointError):#, RuntimeWarning):
		result = 0.
	return result



def expon(x, par):
	return math.exp( (-1.)*float(par[0])*float(x[0]) )

fun = r.TF1("lifetime",expon,0., 160.,2)

def forceDecayInVolume(Momentum, pdfTheta, ct, volume):
	if volume == 1:
		vol = firstVolume
	elif volume == 2:
		vol = secondVolume
	else:
		print "ERROR: select decay volume 1 or 2"
		return 0
	#Direction = Momentum.Vect().Unit()
	#xstart = (-1.)*vol[2] + r.gRandom.Uniform() * (2.*vol[2])
	gamma = Momentum.Gamma()

	fun.SetParameter(0, 1./(gamma*ct))
	
	EndVertex = r.TVector3(0.,0.,vol[1]+1.)
	while EndVertex.Z() >= vol[1]:
		xstart = (-1.)*vol[2] + ((pdfTheta.GetRandom()+math.pi) / (2.*math.pi)) * (2.*vol[2]) 
		Origin = r.TVector3(xstart, 0., vol[0])
		DL = fun.GetRandom(0., 40.)
		EndVertex = r.TVector3(Momentum.Vect().Unit()[0]*DL, Momentum.Vect().Unit()[1]*DL, Momentum.Vect().Unit()[2]*DL) 
		EndVertex = Origin+EndVertex
	#print EndVertex.X(), EndVertex.Z()
	return EndVertex


def makeVtx4AcceptedParaphoton(Momentum, ct, volume):
	if volume == 1:
		vol = firstVolume
	elif volume == 2:
		vol = secondVolume
	else:
		print "ERROR: select decay volume 1 or 2"
		return 0
	gamma = Momentum.Gamma()
	#print Momentum.Px(), Momentum.Py(), Momentum.Pz(), Momentum.E(), ct, gamma*ct
	Direction = Momentum.Vect().Unit()
	costheta = math.fabs(Momentum.CosTheta())
	dxdz = Momentum.Px()/Momentum.Pz()
	dydz = Momentum.Py()/Momentum.Pz()
	Origin = r.TVector3( vol[0]*dxdz, vol[0]*dydz, vol[0] )
	fun.SetParameter(0, 1./(gamma*ct))
	maxlength = 40.*costheta
	DL = fun.GetRandom(0., maxlength)
	EndVertex = r.TVector3(Direction[0]*DL, Direction[1]*DL, Direction[2]*DL)
	EndVertex = Origin + EndVertex
	#print Origin[0], Origin[1], Origin[2]
	#while (EndVertex[0]**2. + EndVertex[1]**2.) > vol[2]:
	#	DL = fun.GetRandom(0., maxlength)
	#	EndVertex = r.TVector3(Direction[0]*DL, Direction[1]*DL, Direction[2]*DL)
	#	EndVertex = Origin + EndVertex
	return EndVertex




def readPDGtable():
    """ Returns R data from PDG in a easy to use format """

    ecm,ratio = [],[]
    with open('rpp2012-hadronicrpp_page1001.dat','r') as f:
        for line in f:
            line = line.split()
            try:
                numEcm = float(line[0])
                numR = float(line[3])
                ecm.append(numEcm)
                ratio.append(numR)
                #print numEcm,numR
            except:
                continue
    #ecm = np.array(ecm)
    #ratio = np.array(ratio)
    return ecm,ratio
    

def interpolatePDGtable(dataEcm,dataR):
    """ Find the best value for R for the given center-of-mass energy """
    fun = interpolate.interp1d(dataEcm,dataR)
    return fun


# sigma(e+e- -> hadrons) / sigma(e+e- -> mu+mu-)
dataEcm,dataR = readPDGtable()
PdgR = interpolatePDGtable(dataEcm,dataR)


#def GeometricAcceptance(gamma, ct, th):
#	if math.fabs(gamma*ct*math.tan(th)) < 2.5 and math.fabs(th) < math.pi/2.:
#		return True
#	return False
def GeometricAcceptance(px, pz, volume):
	if volume == 1:
		vol = firstVolume
	elif volume == 2:
		vol = secondVolume
	else:
		print "ERROR: select decay volume 1 or 2"
		return 0
	#print vol[1]
	#pt = math.sqrt(px*px + py*py)
	# Check theta
    #if (math.fabs((px/pz)*vol[1])<vol[2]) and (pz>0):
    if (math.fabs((px/pz)*vol[0])<vol[2]) and (pz>0):
		return True
	return False


def prodRateFromMesons(mass):
	if pi0Start < mass < pi0Stop:
		chi = 5.41
	elif etaStart < mass < etaStop:
		chi = 0.2314
	elif omegaStart < mass < omegaStop:
		chi = 0.0694
	elif eta1Start < mass < eta1Stop:
		chi = 0.00141
	else:
		print "ERROR: invalid A' mass for meson decay production"
		return 0
	return chi


# Functions for meson decay production
def getTree(mass):
	outfile = r.TFile("out/PythiaData/mesonDecays_%s.root"%mass,"read")
	tree = outfile.Get("mesonDecays")
	return tree


def geomAcceptance(mass, volume):
	if volume == 1:
		tmax = v1ThetaMax
	elif volume == 2:
		tmax = v2ThetaMax
	else:
		print "ERROR: select fiducial volume 1 or 2"
		return -1
	if not os.path.isfile("out/PythiaData/mesonDecays_%s.root"%mass):
		print "ERROR: could not find tree for m(A') = %s"%mass
		return -1
	outfile = r.TFile("out/PythiaData/mesonDecays_%s.root"%mass,"read")
	tree = outfile.Get("mesonDecays")
	#tree = getTree(mass)
	#print tree
	#if not tree:
	#	print "ERROR: could not find tree for m(A') = %s"%mass
	nTot = tree.GetEntriesFast()
	nAcc = 0
	for event in tree:
		if math.fabs(event.A_Theta) < tmax:
			nAcc += 1
	return nAcc / nTot




#SAMPLE USAGE
if __name__ == '__main__':
    ""
    