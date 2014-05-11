from __future__ import division
import sys
import ROOT as r
import math
import numpy as np
import os
import gc
import matplotlib.pyplot as plt
from random import randint
from array import array

from settings import *
from proton_bremsstrahlung import *

upper = "upper"
lower = "lower"

# idea to conserve density of points:
#a = 1.6
#delta = 0.11
#min([i for i in [a + n*delta for n in xrange(100)] if int(i) == i ])


# slope values
alpha0 = -1.*1.5/2.
alpha = -1.*1./2.
beta = -1.*1.2#*2./1.
gamma = 1./0.5#0.5
delta = -1.*1./2.5

# xstart values
a = -3.
b = math.log10(2.e-01)
c = math.log10(4.e-01)
d = math.log10(1.)#0.
e = 1.8#0.6

# ystart values ([upper, lower])
ya = [math.log10(3.e-02), math.log10(2.e-04)]
yb = [math.log10(6.e-05), math.log10(2.e-06)]#[math.log10(6.e-05), math.log10(2.e-06)]
yc = [-7.5, math.log10(9.e-09)]#[math.log10(7.e-08), math.log10(9.e-09)]
yd = [math.log10(1.e-06), math.log10(9.e-08)]#[math.log10(2.e-06), math.log10(9.e-08)]
ye = [-5.8, -8.0]


# new belts (only 2) after seeing the actual shape
# domain: (a, e) = (-3.0, 1.8)
# slopes
sl1 = (-7.5+2.4)/(1.3+3.0)
sl2 = (-8.0+6.6)/(1.3+3.0)
# ystart (upper, lower)
newy1 = [-1.9, -2.9]
newy2 = [-6.1, -7.1]


def eq(n1, n2, toll = 1.e-01):
	if math.fabs(n2-n1) < toll*math.fabs(n2):
		return True
	return False


def lineForTwoPoints(p1, p2, x):
	x1 = p1[0]
	x2 = p2[0]
	y1 = p1[1]
	y2 = p2[1]
	slop = (y2-y1)/(x2-x1)
	return y1 + slop*(x - x1)


def roundToN(x, n=2):
	if x:
		result = round(x, -int(math.floor(math.log10(math.fabs(x)))) + (n - 1))
	else:
		result = 0.
	return result


def beltSegment(pos, ystart, xstart, slope, x):
	if pos == "upper":
		start = ystart[0]
	elif pos == "lower":
		start = ystart[1]
	else:
		print "ERROR: select upper or lower segment!"
	return start - slope*xstart + slope*x



def inBelt((x,y)):
	# old belt segments:
	#if ( (a < x <= b) and (beltSegment(lower, ya, a, alpha, x) < y < beltSegment(upper, ya, a, alpha, x))
	#	or (a < x <= b) and (beltSegment(lower, ya, a, alpha0, x) < y < beltSegment(upper, ya, a, alpha0, x)) ):
	#	return 1
	#elif (b < x <= d) and (beltSegment(lower, yb, b, beta, x) < y < beltSegment(upper, yb, b, beta, x)):
	#	return 2
	#elif (c < x <= d) and (beltSegment(lower, yc, c, gamma, x) < y < beltSegment(upper, yc, c, gamma, x)):
	#	return 3
	#elif (a < x <= c) and (beltSegment(lower, yd, a, delta, x) < y < beltSegment(upper, yd, a, delta, x)):
	#	return 4
	#elif (d < x <= e) and (ye[1] < y < ye[0]):
	#	return 5
	# new belt segments:
	if (a < x < e) and (beltSegment(lower, newy1, a, sl1, x) < y < beltSegment(upper, newy1, a, sl1, x)):
		return True
	# no need for other points in the bottom contour
	#elif (a < x < e) and (beltSegment(lower, newy2, a, sl2, x) < y < beltSegment(upper, newy2, a, sl2, x)):
	#	return True
	return False


def makeSensitivityBelt(existingData, ndivx, ndivy=240, verbose=0):
	points = [(x,y) for x in np.linspace(a, e, ndivx) for y in np.linspace(-9.1, -2.5, ndivy) if inBelt((x,y))]
	data = []
	for i,point in enumerate(points):
		found = False
		mass = roundToN( pow(10.,point[0]) )
		eps  = roundToN( pow(10.,point[1]) )
		for oldDatum in existingData:
			if eq(mass, oldDatum[0]) and eq(eps, oldDatum[1]):
				found = True
				n = oldDatum[2]
				break
		if not found:
			n = roundToN( computeNEvents(mass, eps) )
		logmass = roundToN(point[0])
		logeps = roundToN(point[1])
		datum = ( logmass, logeps, n)
		if verbose:
			print "Point %s of %s: \t log10(mass) %s \t log10(eps) %s \t\t Events: %s"%(i, len(points), logmass, logeps, n)
		data.append(datum)
		gc.collect()
	return data


def makeToySensitivityBelt(ndiv):
	points = [(x,y) for x in np.linspace(a, d, 25) for y in np.linspace(-8., -2.5, 50) if inBelt((x,y))]
	data = []
	for point in points:
		n = float(randint(0, 5))
		datum = (point[0], point[1], n)
		data.append(datum)
	return data


def inTopHalf(p1,p2,datum):
	if datum[1] > lineForTwoPoints(p1,p2,datum[0]):
		return True
	return False


def isInTopContour(datum):
	x = datum[0]
	y = datum[1]
	if ( (a < x <= b) and ( (y > beltSegment(lower, ya, a, alpha, x))
		or (y > beltSegment(lower, ya, a, alpha0, x) ) ) ):
		return True
	elif (b < x <= d) and (y > beltSegment(lower, yb, b, beta, x)):
		return True
	elif ((0.<datum[2]<1000.) and not isInBotContour(datum) and not isInRightContour(datum)):
		print "could not assign datum to contour: ", datum
	return False


def isInBotContour(datum):
	x = datum[0]
	y = datum[1]
	if (c < x <= d) and (y < beltSegment(upper, yc, c, gamma, x)):
		return True
	elif (a < x <= c) and (y < beltSegment(upper, yd, a, delta, x)):
		return True
	return False



def isInRightContour(datum):
	if inBelt((datum[0],datum[1])) == 5:
		return True
	return False


def closest(data,m):
	useful = [dat for dat in data if (m/2. <= dat[2] <= 2.*m)]
	if useful:
		result = min(useful, key=lambda x: math.fabs(x[2]-m))
		return result
	else:
		return False
	
def eq(n1, n2, toll = 1.e-01):
	if math.fabs(n2-n1) < toll*math.fabs(n2):
		return True
	return False


def makeCountours(data,m):
	p1 = (-3., -5.)
	p2 = (1.63, -7.7)
	xAxis = []
	for datum in data:
		xAxis.append(datum[0])
	xAxis = list(set(xAxis))
	topContour = [datum for datum in data if inTopHalf(p1,p2,datum)]
	botContour = [datum for datum in data if not inTopHalf(p1,p2,datum)]
	# Now for every mass value select only the eps value with N closest to m
	topContour = sorted(topContour, key=lambda x: x[0]) #sorted by mass
	botContour = sorted(botContour, key=lambda x: x[0]) #sorted by mass
	top = []
	bot = []
	right = []
	for x in xAxis:
		tempTop = []
		for datum in topContour:
			if datum[0] == x:
				tempTop.append(datum)
		tempBot = []
		for datum in botContour:
			if datum[0] == x:
				tempBot.append(datum)
		if tempTop:
			bestPointTop=closest(tempTop, m)
			if bestPointTop:
				top.append((x, bestPointTop[1], bestPointTop[2]))
		if tempBot:
			bestPointBot=closest(tempBot, m)
			if bestPointBot:
				bot.append((x, bestPointBot[1], bestPointBot[2]))
	top = list(set(top))
	bot = list(set(bot))
	return top, bot



def loadDataFile():
	if not os.path.isfile("out/TextData/sensitivityScan.txt"):
		return 0
	data = []
	with open("out/TextData/sensitivityScan.txt","r") as ifile:
		for line in ifile:
			line = line.split()
			data.append( ( roundToN(float(line[0])),
				roundToN(float(line[1])),
				roundToN(float(line[-1])) ) )
	return data


def convertToLog(data):
	converted = []
	for datum in data:
		converted.append( ( roundToN(math.log10(datum[0])),
			roundToN(math.log10(datum[1])),
			roundToN(datum[2]) ) )
	return converted


def killDuplicates(data):
	comf = data[:]
	for datum in comf:
		copy = data[:]
		copy.remove(datum)
		for check in copy:
			if (eq(datum[0],check[0]) and eq(datum[1], check[1])):
				data.remove(datum)
				break
	return data




if __name__ == '__main__':

	existingData = loadDataFile()
	data = makeSensitivityBelt(existingData, 80, 240, True)
	""" arrive to 0.6 with 30 divisions
	need to set up right contour """
	existingData = convertToLog(existingData)
	data.extend(existingData)
	data = list(set(data))
	data.sort(key=lambda x: x[1])
	data.sort(key=lambda x: x[0])
	top, bottom = makeCountours(data,2.3)
	top.sort(key=lambda x: x[1])
	top.sort(key=lambda x: x[0])
	bottom.sort(key=lambda x: x[1])
	bottom.sort(key= lambda x:-x[0])
	num = len(top) + len(bottom)
	grtot = r.TGraph(num)
	for i in xrange(len(top)):
		grtot.SetPoint(i,top[i][0],top[i][1])
	for k in xrange(len(bottom)):
		grtot.SetPoint(len(top)+k, bottom[k][0],bottom[k][1])
	grtot.SetMarkerStyle(20)
	c1 = r.TCanvas()
	grtot.Draw("alp")
	grtot.GetXaxis().SetTitle(r"log_{10}M_{A}")
	grtot.GetXaxis().SetRangeUser(-3.,2.)
	grtot.GetYaxis().SetTitle(r"log_{10}#varepsilon")
	c1.SetGrid()
	whole = top+bottom
	c2 = r.TCanvas()
	c2.cd()
	grtotlog = r.TGraph(len(whole))
	for i in xrange(len(whole)):
		grtotlog.SetPoint(i, pow(10.,whole[i][0]), pow(10.,whole[i][1]))
	grtotlog.SetMarkerStyle(20)
	grtotlog.Draw()
	grtotlog.GetXaxis().SetTitle(r"M_{A} (GeV)")
	grtotlog.GetYaxis().SetTitle(r"Coupling #varepsilon")
	c2.SetGrid()
	c2.SetLogx()
	c2.SetLogy()

