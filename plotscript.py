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
import copy

global psTimeStep
global vgroup
global cc

psTimeStep = 0.001 #ns
vgroup = 0.218#m/ns
cc = 0.299792458 #m/ns


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
			tempRay = copy.copy(r)
			tempRay.setDetectionTime(tsm)
			tempRay.setDetectionPoint(psm)
			smearedRays.append(tempRay)

	return smearedRays

def dtor(deg):
	return deg*np.pi/180.0

def rtod(rad):
	return rad*180.0/np.pi

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
			z.append(p.getZ()) #in cm
			t.append(r.getDetectionTime())


	ax.scatter(z, t)

	ax.set_xlabel("z along cylinder (m)", fontsize=26)
	ax.set_ylabel("time of detection (ns)", fontsize=26)
	majorLocatorY = MultipleLocator(5)
	minorLocatorY = MultipleLocator(1)
	majorLocatorX = MultipleLocator(.100)
	minorLocatorX = MultipleLocator(.010)
	ax.get_xaxis().set_major_locator(majorLocatorX)
	ax.get_xaxis().set_minor_locator(minorLocatorX)
	ax.get_yaxis().set_major_locator(majorLocatorY)
	ax.get_yaxis().set_minor_locator(minorLocatorY)
	ax.get_xaxis().set_tick_params(labelsize=26, length=20, width=2, which='major')
	ax.get_xaxis().set_tick_params(length=10, width=2, which='minor')
	ax.get_yaxis().set_tick_params(labelsize=26, length=20, width=2, which='major')
	ax.get_yaxis().set_tick_params(length=10, width=2, which='minor')


#plots time vs. z and assigns
#a color based on number of reflections
def plotColorZTime(rays, fig, ax):
	z = []
	t = []
	bounces = []
	for r in rays:
		if(r.isDetected()):
			p = r.getDetectionPoint()
			z.append(p.getZ()) #in cm
			t.append(r.getDetectionTime())
			bounces.append(r.getnRefl())


	maxR = 10
	minR = 0
	cmap = plt.cm.jet
	cmaplist = [cmap(i) for i in range(cmap.N)]
	cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
	bounds = range(minR, maxR)
	norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

	
	ax.scatter(z, t, c=bounces, cmap=cmap, norm=norm, s=30)


	ax2 = fig.add_axes([0.95, 0.1, 0.03, 0.8])
	cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, spacing='proportional', ticks=bounds, boundaries=bounds, format='%1i')
	ax.set_xlabel("z along cylinder (m)", fontsize=26)
	ax.set_ylabel("time of detection (ns)", fontsize=26)
	ax2.set_ylabel("Number of reflections", size=20)
	majorLocatorY = MultipleLocator(5)
	minorLocatorY = MultipleLocator(1)
	majorLocatorX = MultipleLocator(.400)
	minorLocatorX = MultipleLocator(.020)
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
def plotSixColorZTime(rays):
	fig, ax = plt.subplots(nrows=3, ncols=2, figsize=(35, 20))

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
		arrivalTimeHist(colorg[colIndex], fig, ax[colIndex][1])
		#titles
		ax[colIndex][0].set_title("phi = " + str(phis[j][0]), fontsize=25, loc='left')
		ax[colIndex][1].set_title("phi = " + str(phis[j][0]), fontsize=25, loc='left')

	

def plotColorZRotTime(rays, fig, ax):
	z = []
	tp = []
	bounces = []
	for r in rays:
		if(r.isDetected()):
			p = r.getDetectionPoint()
			z.append(p.getZ()) #in cm
			tp.append(r.getDetectionTime() - z[-1]/(cc))
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
	ax.set_xlabel("z along cylinder (m)", fontsize=26)
	ax.set_ylabel("rotated time of detection (ns)", fontsize=26)
	ax2.set_ylabel("Number of reflections", size=20)
	majorLocatorY = MultipleLocator(5)
	minorLocatorY = MultipleLocator(1)
	majorLocatorX = MultipleLocator(.100)
	minorLocatorX = MultipleLocator(.010)
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
	fig, ax = plt.subplots(nrows=3, ncols=2, figsize=(35, 20))

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



def rotatedTimeHist(rays, fig, ax):
	rotatedTime = []
	for r in rays:
		if(r.isDetected()):
			p = r.getDetectionPoint()
			t = r.getDetectionTime()
			z = p.getZ()
			rt = t - (z/cc)
			rotatedTime.append(rt)

	n, bins, patches = ax.hist(rotatedTime, 100,alpha=0.75)
	ax.set_xlabel("rotated time t - z/c (ns)", fontsize=26)
	ax.set_ylabel("detected number of photons", fontsize=26)
	ax.get_xaxis().set_tick_params(labelsize=26, length=20, width=2, which='major')
	ax.get_xaxis().set_tick_params(length=10, width=2, which='minor')
	ax.get_yaxis().set_tick_params(labelsize=26, length=20, width=2, which='major')
	ax.get_yaxis().set_tick_params(length=10, width=2, which='minor')


	
def arrivalTimeHist(rays,fig,ax):
	ts = []
	for r in rays:
		if(r.isDetected()):
			p = r.getDetectionPoint()
			t = r.getDetectionTime()
			z = p.getZ()
			ts.append(t)

	n, bins, patches = ax.hist(ts, 100,alpha=0.75)
	ax.set_xlabel("arrival time (ns)", fontsize=26)
	ax.set_ylabel("detected number of photons", fontsize=26)
	ax.get_xaxis().set_tick_params(labelsize=26, length=20, width=2, which='major')
	ax.get_xaxis().set_tick_params(length=10, width=2, which='minor')
	ax.get_yaxis().set_tick_params(labelsize=26, length=20, width=2, which='major')
	ax.get_yaxis().set_tick_params(length=10, width=2, which='minor')

def lineFitZTime(rays, fig, ax, plot = False):
	z = []
	t = []
	for r in rays:	
		if(r.isDetected()):
			p = r.getDetectionPoint()
			z.append(p.getZ()) #in cm
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
	tr = []
	for r in rays:	
		if(r.isDetected()):
			p = r.getDetectionPoint()
			z.append(p.getZ()) #in cm
			tr.append(r.getDetectionTime() - z[-1]/(cc))

	#fit 1 degree polynomial
	#p = [m, b]
	p, res, _, _, _ = np.polyfit(z, tr, 1, full=True)
	if(plot == True):
		fitz = z
		fitt = [(j*p[0] + p[1]) for j in fitz]
		ax.plot(fitz, fitt, 'g-', linewidth=3, label="dt'/dz = " + str(round(p[0], 4)))
		ax.legend(loc='right', fontsize=25)
		plotColorZRotTime(rays, fig, ax)
		return p[0]
	else:
		return p[0]





#fits a line to photons in z-t space
#finding dt/dz and implementing the geometrical
#constraints to solve for beam angle. Does this for
#each detector row, and weights them to the angle phi of the
#detector
def findBeamAngle(rays, fig, ax):
	#order the rays by column
	colrays = columnOrderedTimeOrdered(eventrays)
	dtdz = []	#[[phicenter, slope], [phicenter, slope]...]
	xhat = 0
	yhat = 0
	for col in colrays:
		true0 = []
		phis = []
		for r in col:
			if(r.getnRefl() == 0):
				true0.append(r)
				phis.append(r.getDetectionPoint().getPhi())


		#fit zero reflection light
		m = lineFitRotZTime(true0, fig, ax, False)

		avgphi=np.mean(phis)
		xhat += -1*np.arctan(m*vgroup)*np.cos(avgphi)
		yhat += -1*np.arctan(m*vgroup)*np.sin(avgphi)


	thI = np.sqrt(xhat**2 + yhat**2)
	phI = np.arccos(xhat/thI)
	print "Beam Angle Theta: " + str(rtod(thI))
	print "Beam angle phi: " + str(rtod(phI))



		


(eventrays, detector) = loadData("data/angled_n3_N5/xyzbeam_phi0_th0_100pcm.p")
#spreadrays = spaceTimeSpread(eventrays, 7e-4, 7e-4, 7e-4, 0.3)
#qerays = qeReduce(spreadrays, 20)

prompt = []
for r in eventrays:
	if(r.isDetected() and r.getnRefl() <= 1):
		prompt.append(r)

fig, ax = plt.subplots(figsize=(20, 12))
plotColorZTime(prompt, fig, ax)
#findBeamAngle(eventrays, fig, ax)
#plt.savefig("testfig.png",bbox_inches='tight')
plt.show()
sys.exit()

	



#plot something for a bunch of files
dirpath = "data/angled_n3_N3/"
filesonly = [f for f in os.listdir(dirpath) if os.path.isfile(os.path.join(dirpath,f))]
for fl in filesonly:
	if(fl[-1] != "p"):
		#not a pickle file
		continue
	else:
		(eventrays, detector) = loadData(dirpath+fl)
		#spreadrays = spaceTimeSpread(eventrays, 7e-4, 7e-4, 7e-4, 0.3)
		#qerays = qeReduce(spreadrays, 20)
		plotSixColorZTime(eventrays)
		plt.savefig("./plots/angled_n3_N3/ColorZ/" + fl[:-2] + ".png", bbox_inches='tight')
		plt.clf()
		plt.close()

sys.exit()








