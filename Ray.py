import Point

import sys
import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



class Ray:
	def __init__(self, p_init, v_init, vgroup, tinit=None):
		self.p = p_init #current point
		self.p_orig = p_init
		self.v = v_init.normalize() #current velocity
		self.v_orig = v_init.normalize()
		self.vgroup = vgroup
		self.track = [p_init]
		self.nRefl = 0 #number of changes in direction (reflections)
		self.pointsSinceLastReflection = 0 #increments when a new point is calculated (used for finding intersections)
		self.detected = False 				#is True when hit photocathode material
		self.absorbed = False
		self.color = np.random.rand(3,1) 	#choose a random trace color
		self.linestyle = '-'				#all rays start as lines
		self.intersectPoints = None		#structure is [[point, region], [point, region], ...]
		self.detectionTime = None		#when the photon is detected
		self.startTime = tinit
		self.detectionPoint = None		#where the photon is is detected
		self.detectionRegion = None		#which region did the detecting
		self.absorptionPoint = None		#which point did it get absorbed
		self.reflectionPoints = []		#appends reflection points


	def isDetected(self):
		return self.detected

	def isAbsorbed(self):
		return self.absorbed

	def setAbsorbed(self, a):
		self.absorbed = a

	def getPos(self):
		return self.p

	def getOrigP(self):
		return self.p_orig

	def getVel(self):
		return self.v

	def getOrigV(self):
		return self.v_orig

	def getnRefl(self):
		return self.nRefl

	def getTrack(self):
		return self.track

	def clearTrack(self):
		self.track = []
		return

	def clearIntersections(self):
		self.intersectPoints = None
		return

	def getDetectionPoint(self):
		return self.detectionPoint

	def getAbsorptionPoint(self):
		return self.absorptionPoint

	def setDetectionPoint(self, p):
		self.detectionPoint = p

	def getDetectionTime(self):
		return self.detectionTime

	def setDetectionTime(self, t):
		self.detectionTime = t

	def getDetectionRegion(self):
		return self.detectionRegion

	def getTrackSize(self):
		return len(self.track)

	def getStartTime(self):
		return self.startTime

	def setStartTime(self, t):
		self.startTime = t #in ns

	def getTravelDist(self):
		travel = 0 #meters
		for i in range(len(self.track) - 1):
			#look at the vector defined by 
			#two neighboring track points
			di = self.track[i]
			df = self.track[i+1]
			d = Point.Point(df.getX() - di.getX(), df.getY() - di.getY(), df.getZ() - di.getZ(), 1)
			#magnitude of this vector is the associated distance
			travel += d.getMag()

		return travel

	def getIntersectionPoints(self):
		return self.intersectPoints

	def setIntersections(self, inters):
		self.intersectPoints = inters

	def shouldFindIntersections(self):
		if(self.pointsSinceLastReflection > 2 and self.intersectPoints == None):
			return True
		else:
			return False

	#check to see if the current position of
	#ray is within a cube of eps length of
	#one of the intersections
	#THEN act on it based on the region
	#of intersection
	def checkIntersection(self, timestep):
		#make sure there are intersection points defined
		if(self.intersectPoints == None):
			return

		possibleInteractionPoints = []
		for sect in self.intersectPoints:
			s = sect[0]	  #the point on the region
			reg = sect[1] #the region containing (or not containing) the point

			#-----selection of interaction criteria------#

			#(1)if the ray is within epsilon sized cube of 
			#the intersection point
			#(2 ispointinregion) if that point is actually inside the region
			#domain, and not just out at infinity
			#(3) there may be multiple regions which satisfy this condition. 
			#is this the first region in time which satisfies the condition
			eps = self.vgroup*timestep
			if(abs(self.p.getX() - s.getX()) < eps and \
				abs(self.p.getY() - s.getY()) < eps and \
				abs(self.p.getZ() - s.getZ()) < eps and \
				reg.isPointInRegion(s)):
				possibleInteractionPoints.append(sect)
			else:
				continue

		#if multiple interaction points are within epsilon
		#and are of a passing region
		#find the one which the ray hits first in time
		interactionPoint = None
		if(len(possibleInteractionPoints) == 0):
			#no interactions were found
			return
		elif(len(possibleInteractionPoints) > 1):
			#multiples found
			([x1, x2], [y1, y2], [z1, z2]) = self.getCurrentLineParameters()
			t = []
			max2 = max([abs(x2),abs(y2),abs(z2)])
			if(max2 < 1e-8):
				#this photon is going "nowhere"
				#remove it arbitrarily, unphysical
				self.absorbed = True
				self.absorptionPoint = self.track[-1]
				return
			else:
				if(max2 == abs(x2)):
					for pr in possibleInteractionPoints:
						x = pr[0].getX()
						t.append((x - x1)/x2)
					
					#minimum t parameter means
					#first intersection point reached
					minT = min(t)

					for pr in possibleInteractionPoints:
						x = pr[0].getX()
						ttemp = (x - x1)/x2
						if(minT == ttemp):
							interactionPoint = pr
				elif(max2 == abs(y2)):
					for pr in possibleInteractionPoints:
						y = pr[0].getY()
						t.append((y - y1)/y2)
					
					#minimum t parameter means
					#first intersection point reached
					minT = min(t)

					for pr in possibleInteractionPoints:
						y = pr[0].getY()
						ttemp = (y - y1)/y2
						if(minT == ttemp):
							interactionPoint = pr
				elif(max2 == abs(z2)):
					for pr in possibleInteractionPoints:
						z = pr[0].getZ()
						t.append((z - z1)/z2)
					
					#minimum t parameter means
					#first intersection point reached
					minT = min(t)

					for pr in possibleInteractionPoints:
						y = pr[0].getZ()
						ttemp = (z - z1)/z2
						if(minT == ttemp):
							interactionPoint = pr

				else:
					print str(max2)
					print self.getCurrentLineParameters()

					print "ROUNDING ERROR (1)"
					sys.exit()
		#only one 
		else:
			interactionPoint = possibleInteractionPoints[0]


		if(interactionPoint == None):
			print "Got an arbitrary rejection"
			self.absorbed = True
			return
		s = interactionPoint[0]
		reg = interactionPoint[1]
		if(reg.getMat() == 'mirror'):
			#reflect
			#get the normal
			n = reg.getNormal()
			#flip the normal if its
			#facing the wrong direction

			if(n.getRho() > 0):
				n = n.scale(-1)

			#the dot product n*v gives 1/2 the fraction
			#of the normal which we want to add to 
			#v in order to "reflect" it.
			cosangle = (self.v*n)/self.v.getMag()
			reflection = n.scale(2*cosangle)
			self.v = self.v - reflection

			#erase the last track point
			#make it the true intersection point with the mirror
			self.track[-1] = s
			self.p = s
			#save data
			self.nRefl += 1
			self.intersectPoints = None
			self.pointsSinceLastReflection = 0
			self.reflectionPoints.append(s)


		elif(reg.getMat() == 'absorber'):
			#forget the ray
			self.absorbed = True
			self.intersectPoints = None
			self.absorptionPoint = s
			
		elif(reg.getMat() == 'cathode'):
			#store detector information

			self.detectionPoint = s
			self.detectionRegion = reg
			self.detected = True
			self.track[-1] = s
			self.p = s
			self.intersectPoints = None


			#calculate the exact time that this photon is detected
			#because right now we may be "eps" away from the point
			travelDist = self.getTravelDist()
			self.detectionTime = travelDist/self.vgroup + self.startTime
		else:
			print "ERROR (1)"
			sys.exit()




	def rotateVel(self, rotMat):
		self.v = self.v.rotate(rotMat)

	def incrementRay(self, stepsize):
		self.p = self.p + self.v.scale(stepsize*self.vgroup)
		self.track.append(self.p)
		self.pointsSinceLastReflection += 1

	def plotRay(self, ax):
		x = [p.getX() for p in self.track]
		y = [p.getY() for p in self.track]
		z = [p.getZ() for p in self.track]
		ax.plot(x,y,z, c=self.color)
		ax.scatter(x[-1], y[-1], z[-1], c='r')
		"""
		if(self.intersectPoints == None):
			return
		for p in self.intersectPoints:
			if(p[1].getShape() != 'square' and p[1].getShape() != 'cylinder'): continue
			ax.scatter(p[0].getX(), p[0].getY(), p[0].getZ(), c='k', marker='d')
		"""

	#returns a parametrization of
	#the line defined by the most recent
	#two points of the track
	#x = ax + (bx - ax)t
	#y = ay + (by - ay)t
	#z = az + (bz - az)t
	def getCurrentLineParameters(self):
		#two most recent points
		pi = self.track[-2]
		pf = self.track[-1]

		#structure is a tuple, 3 components
		#for the three equations above
		#each element of the tuple is a two-list
		#[ai, (bi - ai)]
		return ([pi.getX(), (pf.getX() - pi.getX())], \
				[pi.getY(), (pf.getY() - pi.getY())], \
				[pi.getZ(), (pf.getZ() - pi.getZ())])

