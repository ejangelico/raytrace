import Point
import Ray
import Detector
import Region

import sys
import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
import cPickle as pickle 

global psTimeStep
global c
global groupV

psTimeStep = 0.001 #ns
c = 0.299792458 #m/ns
groupV = 0.218 #m/ns, from eric's paper mean group velocity for 370nm mean wavelength

#propagate a list of rays
def propagateRays(detector, rays, timestep, numberOfSteps):
	for i in range(numberOfSteps):
		if(i % 50 == 0):
			print "On step " + str(i) + " out of " + str(numberOfSteps)
		for r in rays:
			if(r.isDetected() or r.isAbsorbed()):
				continue

			r.incrementRay(timestep)
			#if the intersection points are ready to be calculated
			if(r.shouldFindIntersections()):
				detector.findIntersections(r)

			#check if the rays have reached an intersection point
			#to check if the position is near the intersection
			r.checkIntersection(timestep)


	print "FINISHED PROPAGATING RAYS"




#returns a list of rays which is a conical source
#with an inital position p_init
#and opening angle "opening"
#and number of photons "nrays" distributed around 2pi
#and the axis of the cone is directed with vector "direction"
def generateConeSource(opening, nrays, p_init, direction, vgroup, tinit):
	#get radian interval for nrays
	spacing = 2*np.pi/nrays
	phis = np.arange(0, 2*np.pi, spacing) #in radians

	#calculate z/rho ratio for opening angle
	#assumes opening angle is given in degrees
	ratio = openingAngleRatio(opening)

	zaxis = Point.Point(0,0,1,1)

	#generate a list of rays in a cone 
	#about the origin, then rotate each ray
	rays = []
	for ang in phis:
		v_init = Point.Point(ratio,ang,1,0)
		#rotate this velocity to the specified direction
		directionRotMat = zaxis.getRotationMatrix(direction)
		directionVelocity = v_init.rotate(directionRotMat)
		r = Ray.Ray(p_init, directionVelocity, vgroup, tinit)
		rays.append(r)

	return rays

#generates a set of cones which represent
#nphot per cm traveling in a straight line from
#origin, in the direction of "direction", with velocity
#"particle velocty" and each cone has opening angle of "angle"
def generateMuonEvent(angle,nphot,origin,direction,particleV, trackLength):
	#get the time spacing between cm steps
	#dx = 0.01 #m
	dx = 0.01
	dt = dx/particleV #ns

	#how many distance steps
	nsteps = int(trackLength/dx)

	#make many cones separated in time and space
	rays = []
	print "BUILDING MUON TRACK"
	for i in range(nsteps):
		coneVertex = origin + direction.scale(i*dx)
		coneDirection = direction
		cone = generateConeSource(angle, nphot, coneVertex, coneDirection, groupV, i*dt)
		[rays.append(x) for x in cone]

	return rays




#calculates the ratio of z/rho vector components
#for an opening angle of "theta" in degrees
#will not work for theta in radians
def openingAngleRatio(theta):
	thR = dtor(theta)
	return np.tan(thR)

def dtor(deg):
	return deg*np.pi/180.0

def rtod(rad):
	return rad*180.0/np.pi

#-----Data Analysis----#

def saveData(rays, detector, filename):
	print "saving..."
	smallrays = []
	for r in rays:
		temp = r
		temp.clearIntersections()
		temp.clearTrack()
		smallrays.append(temp)
	pickle.dump([smallrays,detector], open(filename,"wb"))
	print "SAVED DATA to : " + str(filename)
	return

def plotAll(detector, rays):
	print "Drawing detector and event"
	fig = plt.figure(figsize=(11,10))
	ax = fig.gca(projection='3d')
	detector.drawDetector(ax)
	regs = detector.getRegions()
	"""
	for reg in regs:
		if(reg.getShape() == 'square'):
			n = reg.getNormal()
			point = reg.args[1]
			point = Point.Point(0, 0, 0, 0)
			x = [point.getX(), point.getX() + n.getX()/10.]
			y = [point.getY(), point.getY() + n.getY()/10.]
			z = [point.getZ(), point.getZ() + n.getZ()/10.]
			ax.plot(x, y, z, c='r', linewidth=3)
	"""
	[r.plotRay(ax) for r in rays]
	Rad = detector.getR()
	ax.set_xlim([-3*Rad, 3*Rad])
	ax.set_ylim([-3*Rad, 3*Rad])
	ax.set_zlim([0, 0.1 + detector.getL()])
	plt.show()

def getAllReflections(rays):
	ns = []
	for r in rays:
		ns.append(r.getnRefl())

	return ns

def getAllDetectionTimes(rays):
	ts = []
	for r in rays:
		if(r.isDetected()):
			ts.append(r.getDetectionTime())
		else: continue

	return ts









#---------#GEOMETRY
detector = Detector.Detector()
#24'' is 0.6096 m
#96'' is 2.4384 m

#config1 has three rows of 5 detectors
#center phi's: 	3.51075 radians
#				1.41636 radians
#				5.60515 radians
#<number of lappd columns>, <number of LAPPD's per column>, <radius of cylinder>, <mirrored ends>
detector.loadTiltedV2(3, 3, 0.005, 0.005, 0.3048, True)


#detector testing
muonOrigin = Point.Point(0,0,0.2,1)
particleVelocity = c

muonDirection = Point.Point(0,0,1,0).normalize()
cone = generateConeSource(42, 10, muonOrigin, muonDirection, groupV, 0)
propagateRays(detector, cone, psTimeStep, 20000)

plotAll(detector, cone)

sys.exit()




#------#SOURCES

#z-unaligned source
muonOrigin = Point.Point(0,0,0.1,1)
muonTrackLength = 2.43 #m
nphotons = 100 #per cm
particleVelocity = c
tilt = [dtor(4), dtor(7), dtor(15), 0]
r = 1
phis = [dtor(0), dtor(30)] #for each column
for til in tilt:
	if(til == 0):
		muonDirection = Point.Point(0,0,1,1).normalize()
		beam = generateMuonEvent(42, nphotons, muonOrigin, muonDirection, particleVelocity, muonTrackLength)
		propagateRays(detector, beam, psTimeStep, 30000)
		saveData(beam,detector,"data/xyzbeam_phi" + str(0) + "_th" + str(0) + "_100pcm.p")
		beam = 0
	else:
		for ph in phis:
			x = r*np.sin(til)*np.cos(ph)
			y = r*np.sin(til)*np.sin(ph)
			z = r*np.cos(til)
			muonDirection = Point.Point(x,y,z,1).normalize()
			beam = generateMuonEvent(42, nphotons, muonOrigin, muonDirection, particleVelocity, muonTrackLength)
			propagateRays(detector, beam, psTimeStep, 30000)
			saveData(beam,detector,"data/xyzbeam_phi" + str(int(rtod(ph))) + "_th" + str(int(rtod(til))) + "_100pcm.p")
			beam = 0


sys.exit()






