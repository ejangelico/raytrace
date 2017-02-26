import Point
import Region

import sys
import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from itertools import combinations

class Detector:
	def __init__(self):
		self.regions = []
		#variables for cylinder detector only
		self.R = None
		self.L = None
		self.mirrorEnds = False


	#creates a cube detector with
	#bottom face at the origin
	def rectangleDetector(self, xd, yd, zd):
		#testing first with all mirrors
		print "Warning: This function .rectangleDetector is deprecated"
		print "Some standard modules may not work as expected"

		#eight points define a rectangle
		#z = 0 plane
		origin = Point.Point(0,0,0,1)
		z0 = Point.Point(xd/2.0, yd/2.0, 0, 1)
		z1 = Point.Point(-xd/2.0, yd/2.0, 0, 1)
		z2 = Point.Point(xd/2.0, -yd/2.0, 0, 1)
		z3 = Point.Point(-xd/2.0, -yd/2.0, 0, 1)
		#z = zd plane
		zz0 = Point.Point(xd/2.0, yd/2.0, zd, 1)
		zz1 = Point.Point(-xd/2.0, yd/2.0, zd, 1)
		zz2 = Point.Point(xd/2.0, -yd/2.0, zd, 1)
		zz3 = Point.Point(-xd/2.0, -yd/2.0, zd, 1)


		#six faces define a rectangular prism
		faces = [[z0,z1,zz0,zz1], \
				[z0,z2,zz0,zz2], \
				[z2,z3,zz2,zz3], \
		 		[z1,z3,zz1,zz3], \
		 		[zz0,zz1,zz2,zz3]]

		bottomFace = [z0,z1,z2,z3]

		for face in faces:
			self.regions.append(Region.Region('square', 'mirror', face))

		self.regions.append(Region.Region('square', 'cathode', bottomFace))

	#creates a cylindrical detector with
	#its face at the origin
	#and central axis being the z axis
	#"mirrorEnds" is True if you want the end caps of the 
	#cylinder to be mirrors
	def cylindricalDetector(self, R, L, mirrorEnds):
		print "Using a cylindrical detector geometry..."
		self.R = R
		self.L = L
		self.mirrorEnds = mirrorEnds

		cylReg = Region.Region('cylinder', 'absorber', [self.R, self.L])
		self.regions.append(cylReg)

	#defining "saved" detectors 
	#for testing and getting practice
	def loadConfig1(self):
		#create a triangle pattern
		#two of the triangle edges are mirrors
		#the other is photocathode
		print "Initializing detector configuration 1: test configuration..."

		#do the end caps
		if(self.mirrorEnds):
			self.regions.append(Region.Region('circle', 'mirror', [self.R, 0]))
			self.regions.append(Region.Region('circle', 'mirror', [self.R, self.L]))
		else:
			self.regions.append(Region.Region('circle', 'absorber', [self.R, 0]))
			self.regions.append(Region.Region('circle', 'absorber', [self.R, self.L]))


		#mirrors at 0, 120, 240
		lappd = 0.219964
		spacing = 0.01
		z = np.arange(0, self.L - lappd, lappd)
		phi = [self.dtor(0), self.dtor(120), self.dtor(240)]
		mat = 'mirror'
		for a in z:
			for p in phi:
				self.makeCylLAPPD(p, a, mat)

		#alternating mirror/cathode at
		#60, 180, 300
		lappd = 0.219964
		spacing = 0.01
		z = np.arange(0, self.L - lappd, lappd)
		phi = [self.dtor(60), self.dtor(180), self.dtor(300)]
		mats = ['mirror', 'cathode']
		for i in range(len(z)):
			for p in phi:
				if(i % 2 == 0):
					self.makeCylLAPPD(p, z[i], mats[0])
				else:
					self.makeCylLAPPD(p, z[i], mats[1])



	#makes an LAPPD sized photocathode
	#whose corner (closest to origin) is
	#at phi and z and r in cylindrical coordinates

	#the edge of the LAPPD is in contact
	#with the walls of the cylindrical detector
	#i.e. it makes a secant of the cylinder
	def makeCylLAPPD(self, phi, z, mat):

		#basically need 4 points to define the 
		#square region. 
		lappd = 0.219964 #meters length of photocathode

		cornerpoint = Point.Point(self.R, phi, z, 0)
		p1 = Point.Point(self.R, phi, z + lappd, 0) #common edge with corner point

		#find another edge that makes an intersection
		#with the cylinder wall, going in the +phi direction
		phi2 = 2*np.arcsin(lappd/(2*self.R))
		p2 = Point.Point(self.R, phi + phi2, z, 0)
		p3 = Point.Point(self.R, phi + phi2, z + lappd, 0)
		self.regions.append(Region.Region('square', mat, [cornerpoint, p1, p2, p3]))


	#calculate three relevant coverage fractions
	#(1) (area of photocathode) = f1
	#(2) (area of photocathode)/(area of mirror + photocathode) = f2
	#(3) (area of photocathode)/(area of mirror) = f3
	def coverageFractions(self):
		mirrorArea = 0
		pcArea = 0
		for reg in self.regions:
			mat = reg.getMat()
			if(mat == 'mirror'):
				mirrorArea += reg.getArea()
			elif(mat == 'cathode'):
				pcArea += reg.getArea()
			else:
				continue

		f1 = pcArea	#m^2
		f2 = pcArea/(pcArea + mirrorArea)
		f3 = pcArea/mirrorArea

		return (f1, f2, f3)


	
	def drawDetector(self,ax):
		for reg in self.regions:
			reg.drawRegion(ax)

	def findIntersections(self, ray):
		intersections = []
		eps = 1e-6 #micron
		for reg in self.regions:
			sect = reg.findIntersectionPoint(ray)
			if(sect == None):
				continue
			elif(isinstance(sect.x0, complex) or isinstance(sect.x1, complex) or sect.getRho() > (self.R + eps)):
				#point is outside cylindrical detector
				#the equation gives complex intersection point sometimes if outside of cylinder
				continue
			else:
				intersections.append([sect, reg])

		
		ray.setIntersections(intersections)

	def dtor(self,deg):
		return deg*np.pi/180.0

	def rtod(self,rad):
		return rad*180.0/np.pi


