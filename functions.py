import math
import ROOT as r
import numpy as np
from array import array
from scipy import interpolate 
import matplotlib.pyplot as plt

#from settings import *

# useful functions
def energy(p,m):
	""" Compute energy from momentum and mass """
	return math.sqrt(p*p + m*m)

def momentum(E,m):
	""" Compute momentum from energy and mass """
	return math.sqrt(E*E - m*m)

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



#SAMPLE USAGE
if __name__ == '__main__':
    ""
    