import Point
import Ray
import Detector
import Region

import sys
import os
import numpy as np 
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import time
import cPickle as pickle 
from scipy import optimize
from itertools import combinations
import random

global psTimeStep
global speedOfLight

psTimeStep = 0.001 #ns
speedOfLight = 0.218#m/ns


def loadData(filename):
	print "Loading " + filename
	now = time.time()
	f = pickle.load(open(filename, "rb"))
	after = time.time()
	print "Finished loading, took " + str(after-now) + " seconds"
	return (f[0], f[1])

#assumes the last file's detector configuration
#detector is only used in the Plot all function
#so this should be ok
def loadAllData(dirpath):
	#loads all *.p files in the directory given
	filesonly = [f for f in os.listdir(dirpath) if os.path.isfile(os.path.join(dirpath,f))]
	det = None
	bigrays = []
	for fl in filesonly:
		if(fl[-1] != "p"):
			#not a pickle file
			continue
		else:
			(rays, det) = loadData(dirpath+fl)
			[bigrays.append(r) for r in rays]

	return (bigrays, det)

#remove rays that get destroyed by 
#quantum efficiency (constant fraction)
def qeReduce(rays, QE):
	#QE is like "12%" is
	#QE = 12
	reducedRays = []
	for r in rays:
		if(random.randrange(100) < QE):
			reducedRays.append(r)

	return reducedRays

#uses a gaussian distribution spread with sigma
#given by the arguments to change the detected
#position of the rays 
def spaceTimeSpread(rays, xsig, ysig, zsig, tsig):
	smearedRays = []
	for r in rays:
		if(r.isDetected()):
			#get detection coordinates
			t = r.getDetectionTime()
			p = r.getDetectionPoint()
			x = p.getX()
			y = p.getY()
			z = p.getZ()

			#get smeared values
			tsm = np.random.normal(t, tsig)
			xsm = np.random.normal(x, xsig)
			ysm = np.random.normal(y, ysig)
			zsm = np.random.normal(z, zsig)
			psm = Point.Point(xsm, ysm, zsm, 1)

			#make a new set of rays
			tempRay = r
			tempRay.setDetectionTime(tsm)
			tempRay.setDetectionPoint(psm)
			smearedRays.append(tempRay)

	return smearedRays




def plotAll(detector, rays):
	print "Drawing detector and event"
	fig = plt.figure(figsize=(10,8))
	ax = fig.gca(projection='3d')
	detector.drawDetector(ax)
	[r.plotRay(ax) for r in rays]
	plt.show()


def dtor(deg):
	return deg*np.pi/180.0

def rtod(rad):
	return rad*180.0/np.pi

def getAllReflections(rays):
	ns = []
	for r in rays:
		ns.append(r.getnRefl())

	return ns

def getAllDetectionTimes(rays):
	ts = []
	for r in rays:
		if(r.isDetected()):
			if(r.getDetectionTime() < 20):
				ts.append(r.getDetectionTime())
		else: continue

	return ts

def getDetectedFraction(rays):
	tot = len(rays)
	det = 0
	for r in rays:
		if(r.isDetected()):
			det += 1
	return (float(det)/float(tot))

#return an array of detected arrays 
#in time order
def timeOrderedRays(rays):
	detected = []
	for r in rays:
		if(r.isDetected()):
			if(r.getDetectionTime() < 20):
				detected.append(r)

	sortedlist = sorted(detected, key=lambda x: x.getDetectionTime())

	return sortedlist



#orders photons into an array based
#on phi region of detection and time
def regionOrderedTimeOrdered(rays):
	detected = []
	for r in rays:
		if(r.isDetected()):
			detected.append(r)

	sortedRegs = sorted(detected, key=lambda x: x.getDetectionRegion())

	separatedList = []
	templist = []
	currentReg = sortedRegs[0].getDetectionRegion()
	for r in sortedRegs:
		if(r.getDetectionRegion() == currentReg):
			templist.append(r)
		else:
			separatedList.append(templist)
			templist = []
			currentReg = r.getDetectionRegion()
			templist.append(r)

	separatedList.append(templist)

	RorderedTordered = [timeOrderedRays(x) for x in separatedList]
	return RorderedTordered


#orders photons into regions with
#columsn of equal phi
def columnOrderedTimeOrdered(rays):
	detectedList = [r for r in rays if r.isDetected()]
	columnsort = []
	columnsort = sorted(detectedList, key=lambda x: (x.getDetectionRegion()).getPhiBounds())
	columns = []
	temp = []
	currentPhiReg = columnsort[0].getDetectionRegion()
	for colray in columnsort:
		rayreg = colray.getDetectionRegion()
		regphi = rayreg.getPhiBounds()
		if(currentPhiReg.getPhiBounds() == regphi):
			temp.append(colray)
		else:
			columns.append(temp)
			temp = []
			currentPhiReg = rayreg
			temp.append(colray)

	columns.append(temp)
	return columns


def getFirstNDetections(rays, n):
	tordered = timeOrderedRays(rays)
	return tordered[:n]




def plotZTime(rays, fig, ax):
	z = []
	t = []
	for r in rays:
		if(r.isDetected()):
			p = r.getDetectionPoint()
			z.append(p.getZ()*100) #in cm
			t.append(r.getDetectionTime())


	ax.scatter(z, t)

	ax.set_xlabel("z along cylinder (cm)", fontsize=26)
	ax.set_ylabel("time of detection (ns)", fontsize=26)
	majorLocatorY = MultipleLocator(5)
	minorLocatorY = MultipleLocator(1)
	majorLocatorX = MultipleLocator(100)
	minorLocatorX = MultipleLocator(10)
	ax.get_xaxis().set_major_locator(majorLocatorX)
	ax.get_xaxis().set_minor_locator(minorLocatorX)
	ax.get_yaxis().set_major_locator(majorLocatorY)
	ax.get_yaxis().set_minor_locator(minorLocatorY)
	ax.get_xaxis().set_tick_params(labelsize=26, length=20, width=2, which='major')
	ax.get_xaxis().set_tick_params(length=10, width=2, which='minor')
	ax.get_yaxis().set_tick_params(labelsize=26, length=20, width=2, which='major')
	ax.get_yaxis().set_tick_params(length=10, width=2, which='minor')



#makes 6 subplots, 3 columns each with a timing
#distribution and a colorZTime plot
def plotSixColorZTime(rays, filename):
	fig, ax = plt.subplots(nrows=3, ncols=2, figsize=(20, 14))

	#order the columns and plot them consistently
	colorg = columnOrderedTimeOrdered(rays)
	phis = [] #format: [[phi center of column, index in colorg], ...]
	for i in range(len(colorg)):
		reg = colorg[i][0].getDetectionRegion()
		bounds = reg.getPhiBounds()
		phi = 180*(abs((bounds[1] - bounds[0])/2) + min(bounds))/(np.pi)
		phis.append([round(phi, 1), i])

	phis = sorted(phis, key=lambda x: x[1])


	for j in range(len(phis)):
		colIndex = phis[j][1]
		plotColorZTime(colorg[colIndex], fig, ax[colIndex][0])
		#lineFitZTime(colorg[colIndex], fig, ax[colIndex][0], False) #for prompt photons
		arrivalTimeHist(colorg[colIndex], fig, ax[colIndex][1])
		#titles
		ax[colIndex][0].set_title("phi = " + str(phis[j][0]), fontsize=25, loc='left')
		ax[colIndex][1].set_title("phi = " + str(phis[j][0]), fontsize=25, loc='left')

	



#plots time vs. z and assigns
#a color based on number of reflections
def plotColorZTime(rays, fig, ax):
	z = []
	t = []
	bounces = []
	for r in rays:
		if(r.isDetected()):
			p = r.getDetectionPoint()
			z.append(p.getZ()*100) #in cm
			t.append(r.getDetectionTime())
			bounces.append(r.getnRefl())


	maxR = 10
	minR = 0
	cmap = plt.cm.jet
	cmaplist = [cmap(i) for i in range(cmap.N)]
	cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
	bounds = range(minR, maxR)
	norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

	
	ax.scatter(z, t, c=bounces, cmap=cmap, norm=norm, s=100)


	ax2 = fig.add_axes([0.95, 0.1, 0.03, 0.8])
	cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, spacing='proportional', ticks=bounds, boundaries=bounds, format='%1i')
	ax.set_xlabel("z along cylinder (cm)", fontsize=26)
	ax.set_ylabel("time of detection (ns)", fontsize=26)
	ax2.set_ylabel("Number of reflections", size=20)
	majorLocatorY = MultipleLocator(5)
	minorLocatorY = MultipleLocator(1)
	majorLocatorX = MultipleLocator(100)
	minorLocatorX = MultipleLocator(10)
	ax.get_xaxis().set_major_locator(majorLocatorX)
	ax.get_xaxis().set_minor_locator(minorLocatorX)
	ax.get_yaxis().set_major_locator(majorLocatorY)
	ax.get_yaxis().set_minor_locator(minorLocatorY)
	ax.get_xaxis().set_tick_params(labelsize=26, length=20, width=2, which='major')
	ax.get_xaxis().set_tick_params(length=10, width=2, which='minor')
	ax.get_yaxis().set_tick_params(labelsize=26, length=20, width=2, which='major')
	ax.get_yaxis().set_tick_params(length=10, width=2, which='minor')


#makes 6 subplots, 3 columns each with a timing
#distribution and a colorZRotTime plot
def plotSixColorZRotTime(rays, filename):
	fig, ax = plt.subplots(nrows=3, ncols=2, figsize=(25, 12))

	#order the columns and plot them consistently
	colorg = columnOrderedTimeOrdered(rays)
	phis = [] #format: [[phi center of column, index in colorg], ...]
	for i in range(len(colorg)):
		reg = colorg[i][0].getDetectionRegion()
		bounds = reg.getPhiBounds()
		phi = 180*(abs((bounds[1] - bounds[0])/2) + min(bounds))/(np.pi)
		phis.append([round(phi, 1), i])

	phis = sorted(phis, key=lambda x: x[1])


	for j in range(len(phis)):
		colIndex = phis[j][1]
		plotColorZRotTime(colorg[colIndex], fig, ax[colIndex][0])
		#lineFitZTime(colorg[colIndex], fig, ax[colIndex][0], False) #for prompt photons
		rotatedTimeHist(colorg[colIndex], fig, ax[colIndex][1])
		#titles
		ax[colIndex][0].set_title("phi = " + str(phis[j][0]), fontsize=25, loc='left')
		ax[colIndex][1].set_title("phi = " + str(phis[j][0]), fontsize=25, loc='left')


def plotColorZRotTime(rays, fig, ax):
	z = []
	tp = []
	bounces = []
	c = 0.299792458 #m/ns
	for r in rays:
		if(r.isDetected()):
			p = r.getDetectionPoint()
			z.append(p.getZ()*100) #in cm
			tp.append(r.getDetectionTime() - z[-1]/(100*c))
			bounces.append(r.getnRefl())


	maxR = 10
	minR = 0
	cmap = plt.cm.jet
	cmaplist = [cmap(i) for i in range(cmap.N)]
	cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
	bounds = range(minR, maxR)
	norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

	
	ax.scatter(z, tp, c=bounces, cmap=cmap, norm=norm, s=100)


	ax2 = fig.add_axes([0.95, 0.1, 0.03, 0.8])
	cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, spacing='proportional', ticks=bounds, boundaries=bounds, format='%1i')
	ax.set_xlabel("z along cylinder (cm)", fontsize=26)
	ax.set_ylabel("rotated time of detection (ns)", fontsize=26)
	ax2.set_ylabel("Number of reflections", size=20)
	majorLocatorY = MultipleLocator(5)
	minorLocatorY = MultipleLocator(1)
	majorLocatorX = MultipleLocator(100)
	minorLocatorX = MultipleLocator(10)
	ax.get_xaxis().set_major_locator(majorLocatorX)
	ax.get_xaxis().set_minor_locator(minorLocatorX)
	ax.get_yaxis().set_major_locator(majorLocatorY)
	ax.get_yaxis().set_minor_locator(minorLocatorY)
	ax.get_xaxis().set_tick_params(labelsize=26, length=20, width=2, which='major')
	ax.get_xaxis().set_tick_params(length=10, width=2, which='minor')
	ax.get_yaxis().set_tick_params(labelsize=26, length=20, width=2, which='major')
	ax.get_yaxis().set_tick_params(length=10, width=2, which='minor')



def nReflArrivalTime(rays, fig, ax):
	t = []
	n = []
	for r in rays:
		if(r.isDetected()):
			if(r.getDetectionTime() < 20):
				t.append(r.getDetectionTime())
				n.append(r.getnRefl())

	refbins = max(n)

	plt.hist2d(t, n, bins=[40, refbins])
	plt.colorbar()
	ax.set_xlabel("arrival time (ns)", fontsize=26)
	ax.set_ylabel("number of reflections", fontsize=26)
	ax.get_xaxis().set_tick_params(labelsize=26, length=20, width=2, which='major')
	ax.get_xaxis().set_tick_params(length=10, width=2, which='minor')
	ax.get_yaxis().set_tick_params(labelsize=26, length=20, width=2, which='major')
	ax.get_yaxis().set_tick_params(length=10, width=2, which='minor')
	majorLocatorY = MultipleLocator(1)
	majorLocatorX = MultipleLocator(2)
	minorLocatorX = MultipleLocator(0.25)
	ax.get_xaxis().set_major_locator(majorLocatorX)
	ax.get_xaxis().set_minor_locator(minorLocatorX)
	ax.get_yaxis().set_major_locator(majorLocatorY)


def plotPhiZAbsorbed(rays, fig, ax):
	z = []
	phi = []
	for r in rays:
		if(r.isAbsorbed()):
			p = r.getAbsorptionPoint()
			if(p==None):
				continue
			z.append(p.getZ())
			phi.append(p.getPhi()*180/np.pi)

	ax.scatter(phi, z, s=100)

def plotPhiZTime(rays,fig,ax):
	z = []
	t = []
	phi = []
	timeordered = timeOrderedRays(rays)
	for r in timeordered:
			p = r.getDetectionPoint()
			z.append(p.getZ()*100) #in cm
			t.append(r.getDetectionTime())
			phi.append(p.getPhi()*180/np.pi)
	
	maxR = max(t)
	minR = min(t)
	cmap = plt.cm.YlOrRd
	cmaplist = [cmap(i) for i in range(cmap.N)]
	cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
	bounds = np.arange(minR, maxR, (maxR - minR)/10.0)
	norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
	ax.scatter(phi, z, c=t, cmap=cmap, norm=norm, s=100)

	ax2 = fig.add_axes([0.95, 0.1, 0.03, 0.8])
	cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, spacing='proportional', ticks=bounds, boundaries=bounds, format='%3f')
	ax.set_xlabel("phi coordinate (degrees)", fontsize=26)
	ax.set_ylabel("z coordinate (cm)", fontsize=26)
	ax2.set_ylabel("Time of arrival", size=20)
	majorLocatorY = MultipleLocator(100)
	minorLocatorY = MultipleLocator(10)
	majorLocatorX = MultipleLocator(60)
	minorLocatorX = MultipleLocator(10)
	ax.get_xaxis().set_major_locator(majorLocatorX)
	ax.get_xaxis().set_minor_locator(minorLocatorX)
	ax.get_yaxis().set_major_locator(majorLocatorY)
	ax.get_yaxis().set_minor_locator(minorLocatorY)
	ax.get_xaxis().set_tick_params(labelsize=26, length=20, width=2, which='major')
	ax.get_xaxis().set_tick_params(length=10, width=2, which='minor')
	ax.get_yaxis().set_tick_params(labelsize=26, length=20, width=2, which='major')
	ax.get_yaxis().set_tick_params(length=10, width=2, which='minor')


def rotatedTimeHist(rays, fig, ax):
	rotatedTime = []
	c = 0.299792458 #m/ns
	for r in rays:
		if(r.isDetected()):
			p = r.getDetectionPoint()
			t = r.getDetectionTime()
			z = p.getZ()
			rt = t - (z/c)
			rotatedTime.append(rt)

	n, bins, patches = ax.hist(rotatedTime, 30,alpha=0.75)
	ax.set_xlabel("rotated time t - z/c (ns)", fontsize=26)
	ax.set_ylabel("detected number of photons", fontsize=26)
	ax.get_xaxis().set_tick_params(labelsize=26, length=20, width=2, which='major')
	ax.get_xaxis().set_tick_params(length=10, width=2, which='minor')
	ax.get_yaxis().set_tick_params(labelsize=26, length=20, width=2, which='major')
	ax.get_yaxis().set_tick_params(length=10, width=2, which='minor')


#correlates the time of origin (i.e. the
#emittance point of the photon) with the
#z point of detection. 
def plotTimeOriginZ(rays, fig, ax):
	tinit = []
	zdetect = []
	nref = []
	for r in rays:
		if(r.isDetected()):
			p = r.getDetectionPoint()
			zdetect.append(p.getZ()*100) #in cm
			tinit.append(r.getStartTime())
			nref.append(r.getnRefl())


	maxR = 10
	minR = 0
	cmap = plt.cm.jet
	cmaplist = [cmap(i) for i in range(cmap.N)]
	cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
	bounds = range(minR, maxR)
	norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

	
	ax.scatter(tinit, zdetect, c=nref, cmap=cmap, norm=norm, s=100)


	ax2 = fig.add_axes([0.95, 0.1, 0.03, 0.8])
	cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, spacing='proportional', ticks=bounds, boundaries=bounds, format='%1i')
	ax.set_xlabel("time of emittance (ns)", fontsize=26)
	ax.set_ylabel("z of detection (cm)", fontsize=26)
	ax2.set_ylabel("Number of reflections", size=20)
	majorLocatorY = MultipleLocator(100)
	minorLocatorY = MultipleLocator(10)
	majorLocatorX = MultipleLocator(5)
	minorLocatorX = MultipleLocator(1)
	ax.get_xaxis().set_major_locator(majorLocatorX)
	ax.get_xaxis().set_minor_locator(minorLocatorX)
	ax.get_yaxis().set_major_locator(majorLocatorY)
	ax.get_yaxis().set_minor_locator(minorLocatorY)
	ax.get_xaxis().set_tick_params(labelsize=26, length=20, width=2, which='major')
	ax.get_xaxis().set_tick_params(length=10, width=2, which='minor')
	ax.get_yaxis().set_tick_params(labelsize=26, length=20, width=2, which='major')
	ax.get_yaxis().set_tick_params(length=10, width=2, which='minor')


def plotThreeTimeOriginZ(rays):
	fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(15, 15))

	#order the columns and plot them consistently
	colorg = columnOrderedTimeOrdered(rays)
	phis = [] #format: [[phi center of column, index in colorg], ...]
	for i in range(len(colorg)):
		reg = colorg[i][0].getDetectionRegion()
		bounds = reg.getPhiBounds()
		phi = 180*(abs((bounds[1] - bounds[0])/2) + min(bounds))/(np.pi)
		phis.append([round(phi, 1), i])

	phis = sorted(phis, key=lambda x: x[1])


	for j in range(len(phis)):
		colIndex = phis[j][1]
		plotTimeOriginZ(colorg[colIndex], fig, ax[colIndex])
		#titles
		ax[colIndex].set_title("phi = " + str(phis[j][0]), fontsize=25, loc='left')



def rotatedTimeHistColored(rays, fig, ax):
	rotatedTime0 = []
	rotatedTime1 = []
	rotatedTimeUp = []
	c = 0.299792458 #m/ns
	for r in rays:
		if(r.isDetected()):
			p = r.getDetectionPoint()
			t = r.getDetectionTime()
			z = p.getZ()
			rt = t - (z/c)
			if(r.getnRefl() == 0): 
				rotatedTime0.append(rt)
			elif(r.getnRefl() == 1): 
				rotatedTime1.append(rt)
			else:
				rotatedTimeUp.append(rt)

	n, bins, patches = ax.hist(rotatedTime0, 30,alpha=0.75)
	n, bins, patches = ax.hist(rotatedTime1, 30,alpha=0.75)
	n, bins, patches = ax.hist(rotatedTimeUp, 30,alpha=0.75)
	ax.set_xlabel("rotated time t - z/c (ns)", fontsize=26)
	ax.set_ylabel("detected number of photons", fontsize=26)
	ax.get_xaxis().set_tick_params(labelsize=26, length=20, width=2, which='major')
	ax.get_xaxis().set_tick_params(length=10, width=2, which='minor')
	ax.get_yaxis().set_tick_params(labelsize=26, length=20, width=2, which='major')
	ax.get_yaxis().set_tick_params(length=10, width=2, which='minor')

	
def arrivalTimeHist(rays,fig,ax):
	ordered = timeOrderedRays(rays)
	times = getAllDetectionTimes(ordered)
	n, bins, patches = ax.hist(times, 100,alpha=0.75)
	ax.set_xlabel("arrival time (ns)", fontsize=26)
	ax.set_ylabel("detected number of photons", fontsize=26)
	ax.get_xaxis().set_tick_params(labelsize=26, length=20, width=2, which='major')
	ax.get_xaxis().set_tick_params(length=10, width=2, which='minor')
	ax.get_yaxis().set_tick_params(labelsize=26, length=20, width=2, which='major')
	ax.get_yaxis().set_tick_params(length=10, width=2, which='minor')
	majorLocatorX = MultipleLocator(2)
	minorLocatorX = MultipleLocator(0.25)
	ax.get_xaxis().set_major_locator(majorLocatorX)
	ax.get_xaxis().set_minor_locator(minorLocatorX)


#fits each time/z distribution with a line
#for EACH DETECTOR in the set of rays
def lineFitDetectorZTime(rays, fig, ax, plot = False):
	regOrdered = regionOrderedTimeOrdered(rays)
	slopes = [] #slopes of the line fits for each detector
	for region in regOrdered:
		z = []
		t = []
		for r in region:
			if(r.isDetected()):
				p = r.getDetectionPoint()
				z.append(p.getZ()*100) #in cm
				t.append(r.getDetectionTime())

		#fit 1 degree polynomial
		#p = [m, b]
		p, res, _, _, _ = np.polyfit(z, t, 1, full=True)
		if(plot == True):
			fitz = z
			fitt = [(j*p[0] + p[1]) for j in fitz]
			ax.plot(fitz, fitt, '--')
			plotZTime(rays, fig, ax)
		else:
			pass

		slopes.append(p[0])

	return np.average(slopes)

def lineFitZTime(rays, fig, ax, plot = False):
	z = []
	t = []
	for r in rays:	
		if(r.isDetected()):
			p = r.getDetectionPoint()
			z.append(p.getZ()*100) #in cm
			t.append(r.getDetectionTime())

	#fit 1 degree polynomial
	#p = [m, b]
	p, res, _, _, _ = np.polyfit(z, t, 1, full=True)
	if(plot == True):
		fitz = z
		fitt = [(j*p[0] + p[1]) for j in fitz]
		ax.plot(fitz, fitt, 'g-', linewidth=3, label="dt/dz = " + str(round(p[0], 4)))
		ax.legend(loc='right', fontsize=25)
		plotZTime(rays, fig, ax)
		return p[0]
	else:
		return p[0]

def lineFitRotZTime(rays, fig, ax, plot = False):
	z = []
	t = []
	c = 29.9792458 #cm/ns
	for r in rays:	
		if(r.isDetected()):
			p = r.getDetectionPoint()
			z.append(p.getZ()*100) #in cm
			t.append(r.getDetectionTime() - z[-1]/(c))

	#fit 1 degree polynomial
	#p = [m, b]
	p, res, _, _, _ = np.polyfit(z, t, 1, full=True)
	if(plot == True):
		fitz = z
		fitt = [(j*p[0] + p[1]) for j in fitz]
		ax.plot(fitz, fitt, 'g-', linewidth=3, label="dt'/dz = " + str(round(p[0], 4)))
		ax.legend(loc='right', fontsize=25)
		plotColorZRotTime(rays, fig, ax)
		return p[0]
	else:
		return p[0]




#currently guesses the two
#amplitudes of the first distributions
#in the arrival time histogram
#based on the highest count bins
def guessTimeFitParams(n, bins):
	timerange = [min(bins), max(bins)]
	fitTries = int(len(bins)/5)		#number of times to try to fit to
	sigmaRatio = 100 	#rough ratio between A/S
	tries = np.arange(timerange[0], timerange[1], float((timerange[1] - timerange[0])/fitTries))
	a = []
	m = tries
	s = []
	f = []
	tau = []
	l = []
	for t in tries:
		#find closest bin 
		#to the current time
		tbin = None	#index of bins
		for i in range(len(bins) - 1):
			if(t <= bins[i]):
				tbin = i
				break
		if(tbin == None):
			#time not in the range of data
			continue

		a.append(n[tbin])
		s.append(float(n[tbin])/sigmaRatio)

	return (a, m, s)

def fitArbitraryArrivalTime(rays, fig, ax):
	times = getAllDetectionTimes(rays)
	n, bins, patches = ax.hist(times, 140,alpha=0.75)

	fitfunc = lambda p, x: p[0]*np.exp(-1*((x - p[1])**2)/(2*p[2]))
	errfunc = lambda p, x, y: fitfunc(p, x) - y
	a,m,s = guessTimeFitParams(n, bins)
	cbins = [bins[i - 1] + ((bins[i] - bins[i - 1])/2) for i in range(1, len(bins))]
	curves = []
	for i in range(len(a)):
		if(a[i] < 1):
			continue
		p0 = [a[i], m[i], s[i]]
		p1, cov, infodict, mesg, ier = optimize.leastsq(errfunc, p0[:], args=(cbins, n), full_output=True)
		if(p1[1] < 0):
			continue
		if(p1[0] > max(n)):
			continue
		if(p1[1] > max(bins)):
			continue
		else:
			curves.append(p1)

	#throw away duplicates
	#based on "close mean" solutions

	cc = [x for x in combinations(curves, 2)]
	finalFits = []
	for comb in cc:
		fit1 = comb[0]
		fit2 = comb[1]
		#the minimum number of nanoseconds
		#that is allowed for two arrival time
		#distributions is propto the sigma of one of the fits
		minimumSeparation1 = 0.5*np.sqrt(fit1[2])
		minimumSeparation2 = 0.5*np.sqrt(fit2[2])
		meanmeanSeparation = abs(fit1[1] - fit2[1])
		if(meanmeanSeparation < minimumSeparation2 or meanmeanSeparation < minimumSeparation1):
			#remove the later arriving one
			if(fit1[1] > fit2[1]):
				for m in range(len(curves)):
					if(curves[m][1] == fit2[1]):
						curves.pop(m)
						break

			else:
				for m in range(len(curves)):
					if(curves[m][1] == fit1[1]):
						curves.pop(m)
						break
	

	for c in curves:
		#look at the fits
		newT = np.linspace(min(cbins), max(cbins), 100)
		ax.plot(newT, fitfunc(c, newT), linewidth=2)
		ax.plot(c[1] + 2*np.sqrt(c[2]), fitfunc(c, c[1] + 2*np.sqrt(c[2])), 'ro')
		ax.plot(c[1] - 2*np.sqrt(c[2]), fitfunc(c, c[1] - 2*np.sqrt(c[2])), 'ro')
		print "real: " + str(c)
	

def fitAllArrivalTime(rays, fig, ax):
	times = getAllDetectionTimes(rays)
	n, bins, patches = ax.hist(times, 140,alpha=0.75)

	fitfunc = lambda p, x: gausfunc([p[0], p[1], p[2]], x) + gausfunc([p[3], p[4], p[5]], x) + gausfunc([p[6], p[7], p[8]], x) + gausfunc([p[9], p[10], p[11]], x)
	gausfunc = lambda p, x: p[0]*np.exp(-1*((x - p[1])**2)/(2*p[2]))
	errfunc = lambda p, x, y: fitfunc(p, x) - y

	cbins = [bins[i - 1] + ((bins[i] - bins[i - 1])/2) for i in range(1, len(bins))]

	p0 = [650, 0.2, 0.05, 600, 1.7, 0.4, 600, 3.75, 0.4, 400, 8, 4]
	p1, cov, infodict, mesg, ier = optimize.leastsq(errfunc, p0[:], args=(cbins, n), full_output=True)

	newT = np.linspace(min(cbins), max(cbins), 100)
	ax.plot(newT, fitfunc(p1, newT), linewidth=2)
	ax.plot(newT, gausfunc([p1[0], p1[1], p1[2]], newT), linewidth=2, label="mu = " + str(round(p1[1], 3)) + "\ns^2 = " + str(round(p1[2], 3)))
	ax.plot(newT, gausfunc([p1[3], p1[4], p1[5]], newT), linewidth=2, label="mu = " + str(round(p1[4], 3)) + "\ns^2 = " + str(round(p1[5], 3)))
	ax.plot(newT, gausfunc([p1[6], p1[7], p1[8]], newT), linewidth=2, label="mu = " + str(round(p1[7], 3)) + "\ns^2 = " + str(round(p1[8], 3)))
	ax.plot(newT, gausfunc([p1[9], p1[10], p1[11]], newT), linewidth=2, label="mu = " + str(round(p1[10], 3)) + "\ns^2 = " + str(round(p1[11], 3)))
	ax.legend()


#finds the next nearest gaussian looking distribution
#in timing that is no earlier than "earliest" (ns)
#guesses the next clump based on "lastSigma" time away from
#earliest
def fitOneArrivalClump(rays, earliest, lastSigma):
	gausfunc = lambda p, x: p[0]*np.exp(-1*((x - p[1])**2)/(2*p[2]))
	errfunc = lambda p, x, y: gausfunc(p,x) - y

	#time histogram only photons after "earliest"
	validPhotonTimes = []
	for r in rays:
		t = r.getDetectionTime()
		if(t >= earliest):
			validPhotonTimes.append(t)

	n, bins = np.histogram(validPhotonTimes, 100)
	cbins = [bins[i - 1] + ((bins[i] - bins[i - 1])/2) for i in range(1, len(bins))]

	p0 = [max(n), (earliest + 2*lastSigma), lastSigma]
	p1, cov, infodict, mesg, ier = optimize.leastsq(errfunc, p0[:], args=(cbins, n), full_output=True)

	#newT = np.linspace(min(cbins), max(cbins), 100)
	#ax.plot(newT, gausfunc(p1, newT), linewidth=2)
	#plt.show()
	return (p1[1], p1[2]) #returning actual sigma




	


#takes two time ranges:
#range0: the range of time to select what you think is 0 reflections
#range1: ditto, but with 1 reflection
#then takes all photons larger than the upper bound of range1
#and lumps them into >1 reflection
def separateReflections(rays, range0, range1):
	refl0 = []
	refl1 = []
	reflUp = []
	timeord = timeOrderedRays(rays)
	for r in timeord:
		artime = r.getDetectionTime()
		if(artime <= range0[1] and artime >= range0[0]):
			refl0.append(r)
		elif(artime <= range1[1] and artime >= range1[0]):
			refl1.append(r)
		elif(artime > range1[1]):
			reflUp.append(r)
		else:
			#just throw it away
			continue

	return (refl0, refl1, reflUp)

#takes time cut values from the arrival time fit
#(must be hard coded in)
#and separates rays by "guessed" reflection
#then uses simulated information to see how good
#that guess is

#this function should be improved to allow for
#arbitrary number of reflections, so 3 - 8 instead of
#ref0, ref1, refUp
def plotSeparationEfficiency(rays, fig, ax):
	sigma0 = np.sqrt(0.038)	#ns
	mean0 = 0.363	#ns
	sigma1 = np.sqrt(0.892)	#ns
	mean1 = 1.558	#ns

	total0refl = 0
	total1refl = 0
	totaluprefl = 0
	for r in rays:
		if(r.isDetected()):
			nr = r.getnRefl()
			if(nr == 0):
				total0refl += 1
			elif(nr == 1):
				total1refl += 1
			else:
				totaluprefl += 1

	howmanysigma = np.arange(0.1, 6, 0.1)
	eff0 = []
	eff1 = []
	effUp = []
	max0 = [0,0]
	max1 = [0,0]
	for s in howmanysigma:
		range0 = [mean0 - s*sigma0, mean0 + s*sigma0]
		range1 = [mean1 - s*sigma1, mean1 + s*sigma1]
		(ref0, ref1, refup) = separateReflections(rays, range0, range1)
		G0 = len(ref0)
		G1 = len(ref1)
		GUp = len(refup)
		g0 = 0
		g1 = 0
		gUp = 0
		for _ in ref0:
			if(_.getnRefl() == 0):
				g0 += 1
		for _ in ref1:
			if(_.getnRefl() == 1):
				g1 += 1
		for _ in refup:
			if(_.getnRefl() > 1):
				gUp += 1

		selectEff0 = 2*g0/float(G0 + g0)
		selectEff1 = 2*g1/float(G1 + g1)
		selectEffUp = 2*gUp/float(GUp + gUp)
		globalEff0 = 1 - abs(G0 - total0refl)/float(G0 + total0refl)
		globalEff1 = 1 - abs(G1 - total1refl)/float(G1 + total1refl)
		globalEffUp = 1 - abs(GUp - totaluprefl)/float(GUp + totaluprefl)
		eff0.append(selectEff0*globalEff0)
		eff1.append(selectEff1*globalEff1)
		effUp.append(selectEffUp*globalEffUp)
		if(max0[0] < selectEff0*globalEff0):
			max0[0] = selectEff0*globalEff0
			max0[1] = s
		if(max1[0] < selectEff1*globalEff1):
			max1[0] = selectEff1*globalEff1
			max1[1] = s


	ax.plot(howmanysigma, eff0, linewidth=2, label="0 refl efficiencymax = " + str(max0[0]) + " at sig = " + str(max0[1]))
	ax.plot(howmanysigma, eff1, linewidth=2, label="1 refl efficiency, max = " + str(max1[0]) + " at sig = " + str(max1[1]))
	ax.plot(howmanysigma, effUp, linewidth=2, label="More than 1 refl efficiency")
	ax.legend(loc='lower right')
	ax.set_title("Reflection Selection efficiency: (1 - |G-N|/g-n)(2g/G+g)", fontsize=26)
	ax.set_xlabel("sigma selection range (ns)", fontsize=26)
	ax.set_ylabel("efficiency", fontsize=26)
	ax.get_xaxis().set_tick_params(labelsize=26, length=20, width=2, which='major')
	ax.get_xaxis().set_tick_params(length=10, width=2, which='minor')
	ax.get_yaxis().set_tick_params(labelsize=26, length=20, width=2, which='major')
	ax.get_yaxis().set_tick_params(length=10, width=2, which='minor')
	majorLocatorX = MultipleLocator(1)
	minorLocatorX = MultipleLocator(0.25)
	majorLocatorY = MultipleLocator(0.1)
	ax.get_xaxis().set_major_locator(majorLocatorX)
	ax.get_xaxis().set_minor_locator(minorLocatorX)
	ax.get_yaxis().set_major_locator(majorLocatorY)
	ax.grid(b=True, which='major', color='k', linestyle='--', alpha=0.75)


#actually selects 0, 1, and more reflection
#photons based on pre-calculated/simulated 
#selection efficiencies
#returns ref0, ref1, refUp
def selectReflectionPhotons(rays):
	sigma0 = np.sqrt(0.038)	#ns
	mean0 = 0.363	#ns
	sigma1 = np.sqrt(0.892)	#ns
	mean1 = 1.558	#ns

	#calculated/simulated from 
	#PlotSelectionEfficiency
	ref0Sigma = 5.9
	ref1Sigma = 0.7
	range0 = [mean0 - ref0Sigma*sigma0, mean0 + ref0Sigma*sigma0]
	range1 = [mean1 - ref1Sigma*sigma1, mean1 + ref1Sigma*sigma1]
	#mind you, this gives overlap between the sets
	#ref0 and ref1
	return separateReflections(rays, range0, range1)


#expects column separated rays already
#for small angle events, the first few photons
#in the lowest z-detector are almost 100% prompt
#photons. 
#Algorithm:
#grab the clump of photons in the first z-detector
#in the next over z-detector, grab the nearest clump of
#photons that are all arriving later than the latest photon
#from the previous detector
def findFirstResponders(colrays):
	promptCols = []
	for col in colrays:
		regord = regionOrderedTimeOrdered(col)
		zord = sorted(regord, key=lambda x: x[0].getDetectionRegion().getZBounds())
		lastLatestArrival = -2 #ns
		lastSigma = 0.1 #ns
		prompts = []
		for zd in zord[:3]:
			#fig, ax = plt.subplots(figsize=(15, 10))
			mean, sigma = fitOneArrivalClump(zd, lastLatestArrival, lastSigma)
			#save all photons within n*sigma of this clump fit
			saveSigma = 1
			maxtime = -10
			for r in zd:
				t = r.getDetectionTime()
				saveRange = [mean - saveSigma*np.sqrt(sigma), mean + saveSigma*np.sqrt(sigma)]
				if(t <= saveRange[1] and t >= saveRange[0]):
					prompts.append(r)
					if(t > maxtime):
						maxtime = t

			lastLatestArrival = maxtime
			lastSigma = sigma

		promptCols.append(prompts)

	plotsixrays = [item for sublist in promptCols for item in sublist]
	plotSixColorZTime(plotsixrays, "prompts")
	#plt.show()
	plt.savefig("plots/beamrotate/ColorZ/6plot/find_prompts_phi80_th-4.png", bbox_inches='tight')


			




#fits a line to photons in z-t space
#finding dt/dz and implementing the geometrical
#constraints to solve for beam angle. Does this for
#each detector row, and weights them to the angle phi of the
#detector
def findBeamAngle(rays, fig, ax):
	#order the rays by column
	colrays = columnOrderedTimeOrdered(eventrays)
	dtdz = []	#[[phicenter, slope], [phicenter, slope]...]
	c = 29.9792458 #cm/ns
	for col in colrays:
		#make reflection cuts
		#ref0, ref1, refup = selectReflectionPhotons(col)
		true0 = [] 
		for r in col:
			if(r.getnRefl() == 0):
				true0.append(r)

		#fit zero reflection light
		m = lineFitRotZTime(true0, fig, ax, True)
		print "Slope = " + str(m)
		print "implies an angle of " + str(rtod(np.arctan(m*speedOfLight)))
		plt.show()

		

#need to use pen and paper to work out geometry
#but this solves for a beam incident angle 
#given detector geometry and cerenkov angle
#dtdz= [[centerphi, dtdz], [centerphi, dtdz], ... ] in cm/ns
#ch in degrees
def calculateIncAngle(dtdz, ch):
	chrad = dtor(ch)
	betac = 29.97 	#cm/ns (beta*c = c)
	vg = speedOfLight*100 #group velocity of the light in the water


	#numerically solve a serious
	#transcendental function of the incident angle
	#calling incident angle, Ti for "theta_i"
	listTi = []
	for row in dtdz:
		sl = row[1]
		func = lambda ti: ((1.0/betac - (1.0/vg)*np.sin(ti)*(1.0/np.sin(ti + chrad)))/(np.cos(ti) - np.sin(ti)*(1.0/np.tan(chrad + ti))) - sl)**2

		#solve for minimum
		#in between the physical divergences, 0 and 90 degrees
		sol = optimize.minimize_scalar(func, bounds=(-np.pi/2.0 + 1e-6, np.pi/2.0 - 1e-6), method='bounded')
		listTi.append([row[0], sol.x]) #in radians
		print rtod(sol.x)

	#based on the three angles, one for each
	#row of detectors, find the phi around the 
	#cylindrical axis
	resPhi, resTh = incidentPhi(listTi)
	print rtod(resPhi)
	print rtod(resTh)

#based on the three rows of detectors
#calculate the phi about the cylindrical
#axis, using planar incident elevated angles
#as weights in polar coordinates 
def incidentPhi(ti):
	xyzvectors = []
	for row in ti:
		R = row[1] #in "radians"
		phi = row[0] #in radians
		x = R*np.sin(phi)
		y = R*np.cos(phi)
		xyzvectors.append([x,y])

	resultant = [0,0]
	for v in xyzvectors:
		resultant[0] += v[0]
		resultant[1] += v[1]

	rphi = np.arctan(resultant[1]/resultant[0])
	rtheta = np.sqrt(resultant[1]**2 + resultant[0]**2)
	return (rphi, rtheta)


#currently spits out the position and 
#detection times of the first three photons,
#one from each column's first Z detector
def findFirstFewPhotons(rays, phitrue, thtrue):
	#order the columns and plot them consistently
	colorg = columnOrderedTimeOrdered(rays)
	ps = []
	ts = []
	tinits = []
	for col in colorg:
		regord = regionOrderedTimeOrdered(col)
		zord = sorted(regord, key=lambda x: x[0].getDetectionRegion().getZBounds())
		ztord = timeOrderedRays(zord[0])#look at first z-detector
		firstPhoton = ztord[0]
		t = firstPhoton.getDetectionTime()
		p = firstPhoton.getDetectionPoint()
		ps.append(p)
		ts.append(t)
		tinits.append(firstPhoton.getStartTime())

	#solve a simple linear equation to get 
	#a rough idea of the incoming particle direction
	b = []
	a = []
	vg = speedOfLight
	bc = 0.299792458 #m/ns
	for i in range(len(ps)):
		b.append(ts[i]*vg*vg)
		a.append([ps[i].getX(), ps[i].getY(), ps[i].getZ()])

	a = np.array(a)
	b = np.array(b)
	u = np.linalg.solve(a, b)
	su = Point.Point(u[0], u[1], u[2], 1).normalize()
	print "Volocity in (x, y, z) = (" + str(su.getX()) + ", " + str(su.getY()) + ", " + str(su.getZ()) + ")"

	#true value for comparison
	til = dtor(thtrue)
	ph = dtor(phitrue)
	x = 1*np.sin(til)*np.cos(ph)
	y = 1*np.sin(til)*np.sin(ph)
	z = 1*np.cos(til)
	mu = Point.Point(x,y,z,1).normalize()
	print "phi" + str(phitrue) + " and th" +str(thtrue) + " from raytrace = (" + str(mu.getX()) + ", " + str(mu.getY()) + ", " + str(mu.getZ()) + ")"

	#error, just the angle between
	#these vectors
	print "Dot product of reconstructed and true: " + str(rtod(np.arccos(mu*su)))
	#were these photons really originating from the same point?
	print "The emittance times: " 
	for i in range(3):
		print str(tinits[i])

	diffs = [abs(x[1] - x[0]) for x in combinations(tinits, 2)]
	print "Max time difference is " + str(max(diffs))
	print "Gives max distance difference (*bc) = " + str(max(diffs)*bc)








#(allrays, detector) = loadAllData("data/beamrotate/")
(eventrays, detector) = loadData("data/beamrotate/xyzbeam_phi81_th7_100pcm.p")
spreadrays = spaceTimeSpread(eventrays, 7e-4, 7e-4, 7e-4, 0.3)
qerays = qeReduce(spreadrays, 20)
#colrays = columnOrderedTimeOrdered(eventrays)
fig, ax = plt.subplots(figsize=(25, 18))
#plotThreeTimeOriginZ(spreadrays)
#plotSixColorZTime(eventrays, "tse")
#plotSixColorZRotTime(eventrays, "test")
findBeamAngle(qerays, fig, ax)
#prompt = findFirstResponders(colrays)
#findBeamAngle(eventrays, detector, cerenkovAngle, fig, ax)
#plt.savefig("testfig.png",bbox_inches='tight')
plt.show()
sys.exit()

#plot something for a bunch of files
dirpath = "data/beamrotate/"
filesonly = [f for f in os.listdir(dirpath) if os.path.isfile(os.path.join(dirpath,f))]
for fl in filesonly:
	if(fl[-1] != "p"):
		#not a pickle file
		continue
	else:
		(eventrays, detector) = loadData(dirpath+fl)
		spreadrays = spaceTimeSpread(eventrays, 7e-4, 7e-4, 7e-4, 0.3)
		qerays = qeReduce(spreadrays, 20)
		plotSixColorZTime(qerays, fl[:-2])
		plt.savefig("./plots/beamrotate/ColorZ/6plot/qe20_" + fl[:-2] + ".png", bbox_inches='tight')
		plt.clf()
		plt.close()

sys.exit()











	
#------#Calibration procedure for new large data sets
(allrays, detector) = loadAllData("data/beamrotate/")
fig, ax = plt.subplots(figsize=(12, 8))

#check where the 0 and 1 refl distributions peak
#nReflArrivalTime(allrays, fig, ax)

#visually inspect the fit of the timing distribution
#adjust if necessary
#fitAllArrivalTime(allrays, fig, ax)
#hardcode sigmas/means into plotSeparationEfficiency


#plot this and make sure that you are correctly
#identifying first reflection/second reflection time
#distributions
plotSeparationEfficiency(allrays, fig, ax)
#hardcode the maximum sigma into selectReflectionPhotons

#then you're done

plt.show()
sys.exit()
