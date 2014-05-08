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

# slope values
#alpha = -1.*1.5/2.
alpha = -1.*1./2.
beta = -1.*2./1.
gamma = 1./0.5
delta = -1.*1./2.5

# xstart values
a = -3.
b = math.log10(2.e-01)
c = math.log10(4.e-01)
d = math.log10(1.)#0.

# ystart values ([upper, lower])
ya = [math.log10(3.e-02), math.log10(2.e-04)]
yb = [math.log10(6.e-05), math.log10(2.e-06)]
yc = [math.log10(7.e-08), math.log10(9.e-09)]
yd = [math.log10(2.e-06), math.log10(9.e-08)]


def eq(n1, n2, toll = 1.e-01):
	if math.fabs(n2-n1) < toll*math.fabs(n2):
		return True
	return False




def beltSegment(pos, ystart, xstart, slope, x):
	if pos == "upper":
		start = ystart[0]
	elif pos == "lower":
		start = ystart[1]
	else:
		print "ERROR: select upper or lower segment!"
	return start - slope*xstart + slope*x


def inBelt((x,y)):
	if (a < x <= b) and (beltSegment(lower, ya, a, alpha, x) < y < beltSegment(upper, ya, a, alpha, x)):
		return 1
	elif (b < x <= d) and (beltSegment(lower, yb, b, beta, x) < y < beltSegment(upper, yb, b, beta, x)):
		return 2
	elif (c < x <= d) and (beltSegment(lower, yc, c, gamma, x) < y < beltSegment(upper, yc, c, gamma, x)):
		return 3
	elif (a < x <= c) and (beltSegment(lower, yd, a, delta, x) < y < beltSegment(upper, yd, a, delta, x)):
		return 4
	return 0


def makeSensitivityBelt(existingData, ndiv, verbose=0):
	points = [(x,y) for x in np.linspace(a, d, ndiv) for y in np.linspace(-8., -2.5, ndiv*2) if inBelt((x,y))]
	#listx = []
	#listy = []
	#for point in points:
	#	listx.append(point[0])
	#	listy.append(point[1])
	#plt.ion()
	#plt.scatter(listx, listy)
	data = []
	for i,point in enumerate(points):
		found = False
		mass = pow(10.,point[0])
		eps = pow(10.,point[1])
		for oldDatum in existingData:
			if eq(mass, oldDatum[0]) and eq(eps, oldDatum[1]):
			#if eq(point[0], oldDatum[0]) and eq(point[1], oldDatum[1]):
				found = True
				n = oldDatum[2]
				break
		if not found:
			n = computeNEvents(mass, eps)
		datum = (point[0], point[1], n)
		if verbose:
			print "Point %s of %s: \t log10Mass %s \t log10Eps %s \t\t Events: %s"%(i, len(points), point[0], point[1], n)
		data.append(datum)
		gc.collect()
	return data


def makeToySensitivityBelt(ndiv):
	points = [(x,y) for x in np.linspace(a, d, ndiv) for y in np.linspace(-8., -2.5, ndiv*2) if inBelt((x,y))]
	data = []
	for point in points:
		n = float(randint(0, 5))
		datum = (point[0], point[1], n)
		data.append(datum)
	return data


def isInTopContour(datum):
	segment = inBelt((datum[0],datum[1]))
	if segment == 1 or segment == 2:
		return True
	elif segment == 3 or segment == 4:
		return False
	else:
		print "could not assign datum to contour: ", datum
	return False


def closest(data,m):
	#if not data:
	#	print "Error"
	useful = [dat for dat in data if (dat[2]>=m/2.)]
	if useful:
		result = min(useful, key=lambda x: math.fabs(x[2]-m))
		print result
		return result
	else:
		return False
	
def eq(n1, n2, toll = 1.e-01):
	if math.fabs(n2-n1) < toll*math.fabs(n2):
		return True
	return False


def makeCountours(data,m):
	xAxis = []
	for datum in data:
		xAxis.append(datum[0])
	xAxis = list(set(xAxis))
	topContour = [datum for datum in data if isInTopContour(datum)]
	botContour = [datum for datum in data if not isInTopContour(datum)]
	# Now for every mass value select only the eps value with N closest to m
	topContour = sorted(topContour, key=lambda x: x[0]) #sorted by mass
	botContour = sorted(botContour, key=lambda x: x[0]) #sorted by mass
	top = []
	bot = []

	for x in xAxis:
		tempTop = []
		for datum in topContour:
			if datum[0] == x:
				tempTop.append(datum)
		tempBot = []
		for datum in botContour:
			if datum[0] == x:
				tempBot.append(datum)
		#if (not tempBot) or (not tempTop):
		#	print "error at x= ",x
		if tempTop:
			bestPointTop=closest(tempTop, m)
			if bestPointTop:
				top.append((x, bestPointTop[1]))
		if tempBot:
			bestPointBot=closest(tempBot, m)
			if bestPointBot:
				bot.append((x, bestPointBot[1]))
	return top, bot



def loadDataFile():
	if not os.path.isfile("out/TextData/sensitivityScan.txt"):
		return 0
	data = []
	with open("out/TextData/sensitivityScan.txt","r") as ifile:
		for line in ifile:
			line = line.split()
			data.append((float(line[0]),float(line[1]),float(line[-1])))
	return data


def convertToLog(data):
	converted = []
	for datum in data:
		converted.append((math.log10(datum[0]), math.log10(datum[1]), datum[2]))
	return converted


def killDuplicates(data):
	comf = [myitem(el) for       ] data[:]
	for datum in comf:
		copy = data[:]
		copy.remove(datum)
		for check in copy:
			if (eq(datum[0],check[0]) and eq(datum[1], check[1])):
				data.remove(datum)
				break
	return data




if __name__ == '__main__':

	#print alpha, beta, gamma, delta
	#print a, b, c, d
	existingData = loadDataFile()
	data = makeSensitivityBelt(existingData, 25, True)
	existingData = convertToLog(existingData)
	data.extend(existingData)
	#print len(data)
	#data = killDuplicates(data)
	#print len(data)
	data = list(set(data))
	data.sort(key=lambda x: x[1])
	data.sort(key=lambda x: x[0])
	top, bottom = makeCountours(data,2.3)
	#sys.exit(0)
	top = killDuplicates(top)
	bottom = killDuplicates(bottom)
	xt = []
	yt = []
	xb = []
	yb = []
	top.sort(key=lambda x: x[1])
	top.sort(key=lambda x: x[0])
	bottom.sort(key=lambda x: x[1])
	bottom.sort(key=lambda x: x[0])
	for point in top:
		xt.append(point[0])
		yt.append(point[1])
	for point in bottom:
		xb.append(point[0])
		yb.append(point[1])
	plt.ion()
	plt.scatter(xt,yt)
	plt.scatter(xb,yb)
	mg = r.TMultiGraph()
	tgr = r.TGraph(len(top),array('f',xt),array('f',yt))
	bgr = r.TGraph(len(bottom),array('f',xb),array('f',yb))
	tgr.SetMarkerStyle(8)
	bgr.SetMarkerStyle(8)
	mg.Add(tgr)
	mg.Add(bgr)
	mg.Draw("acp")
	mg.GetXaxis().SetTitle(r"log_{10}M_{A}")
	mg.GetYaxis().SetTitle(r"log_{10}#varepsilon")