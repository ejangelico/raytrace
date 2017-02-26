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

def generateIsotropicSource(nrays, p_init, vgroup, tinit):
	bomb = []
	for n in range(nrays):
		phi = np.random.uniform(0, 2*np.pi)
		costheta = np.random.uniform(-1,1)
		theta = np.arccos(costheta)
		r = 1
		x = r*np.sin(theta)*np.cos(phi)
		y = r*np.sin(theta)*np.sin(phi)
		z = r*np.cos(theta)
		direction = Point.Point(x, y, z, 1).normalize()
		bomb.append(Ray.Ray(p_init, direction, vgroup, tinit))

	return bomb



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
		smallrays.append(temp)
	pickle.dump([smallrays,detector], open(filename,"wb"))
	print "SAVED DATA to : " + str(filename)
	return

def plotAll(detector, rays):
	print "Drawing detector and event"
	fig = plt.figure(figsize=(15,15))
	ax = fig.gca(projection='3d')
	ax.set_xlim([-1.5, 1.5])
	ax.set_ylim([-1.5, 1.5])
	ax.set_zlim([-1.5, 1.5])
	detector.drawDetector(ax)
	plt.axis('off')
	regs = detector.getRegions()
	[r.plotRay(ax) for r in rays]
	plt.show()
	#plt.savefig("detectorPics/cube_vertex45.png", bbox_inches='tight')
	

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
detector.mirrorRect(1, 1, 4.188)
#detector.mirrorSphere(1)
#detector.mirrorCylinder(0.8, 2.07)

#detector testing
muonOrigin = Point.Point(0.3,0.3,0.3,1)
direction = Point.Point(-0.3,0,1,1)
particleVelocity = c

cone = generateMuonEvent(42, 30,muonOrigin,direction,particleVelocity,0.3)

propagateRays(detector, cone, psTimeStep, 100000)

#plotAll(detector, cone)
saveData(cone, detector, "data/cube_muon_100ns.p")

sys.exit()







#vertex source
muonDirection = Point.Point(0,0,1,1).normalize()
cone = generateConeSource(42, 50, muonOrigin, muonDirection, groupV, 0)
muonDirection = Point.Point(0,1,1,1).normalize()
cone2 = generateConeSource(42, 50, muonOrigin, muonDirection, groupV, 0)
for c in cone2:
	cone.append(c)




