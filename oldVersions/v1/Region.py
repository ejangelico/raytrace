import Point

import sys
import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from itertools import combinations
import cPickle as pickle

class Region:
	def __init__(self, shape, mat, args):
		#shape is a string which defines the 
		#shape of the desired patch

		#args is a list of arguments
		#to pass to the shape building functions
		self.shape = shape


		#different for each shape:
		#square [point1, point2, point3, point4] defining the square
		#circle [radius, zposition along axis] (only cirles on axis)
		#cylinder [radius, length] (aligned on the axis with face on origin)
		self.args = args
		self.mat = mat
		self.color = self.getColor()


	def __str__(self):
		st = "Shape: " + str(self.shape) + "\n"
		st += "Mat: " + str(self.mat) + "\n"
		st += "Args: " 
		st += str(self.args) + "\n"
		return st

	def getColor(self):
		if (self.mat == 'mirror'):
			return 'c'
		elif (self.mat == 'absorber'):
			return 'k'
		elif (self.mat == 'cathode'):
			return 'r'

		else:
			print "Material is incorrect! Check detector definition"
			return 'k'

	def getMat(self):
		return self.mat

	def setMat(self, mat):
		self.mat = mat

	def getShape(self):
		return self.shape

	def drawRegion(self, ax):
		if(self.shape == 'square'):
			self.drawSquare(ax)
		elif(self.shape == 'circle'):
			self.drawCircle(ax)
		elif(self.shape == 'cylinder'):
			self.drawCylinder(ax)


	def drawSquare(self, ax):
		#draw all lines connecting the four 
		#points in the argument, which
		#represent the four corners

		lines = [x for x in combinations(self.args, 2)]
		for line in lines:
			x = [p.getX() for p in line]
			y = [p.getY() for p in line]
			z = [p.getZ() for p in line]
			ax.plot(x,y,z, c=self.color)

	def drawCircle(self, ax):
		#arguments for a circle are
		#args = [radius, z-position]
		zshift = self.args[1]
		r = self.args[0]
		plotpoints = []
		npoints = 30.0
		thstep = np.pi/npoints #radians
		angs = np.arange(0, 2*np.pi, thstep)
		for a in angs:
			plotpoints.append(Point.Point(r, a, zshift, 0))

		x = [p.getX() for p in plotpoints]
		y = [p.getY() for p in plotpoints]
		z = [p.getZ() for p in plotpoints]

		ax.plot(x, y, z, color=self.getColor())

	def drawCylinder(self, ax):
		nlines = 8.0
		length = self.args[1]
		dphi = 2*np.pi/nlines
		r = self.args[0]
		plotpoints = []
		angs = np.arange(0, 2*np.pi, dphi)
		for a in angs:
			plotpoints.append(Point.Point(r, a, 0, 0))
			plotpoints.append(Point.Point(r, a, length, 0))
			x = [p.getX() for p in plotpoints]
			y = [p.getY() for p in plotpoints]
			z = [p.getZ() for p in plotpoints]
			ax.plot(x, y, z, color=self.getColor(), linestyle='-')
			plotpoints = []

	def getArea(self):
		if(self.shape == 'square'):
			lines = [x for x in combinations(self.args, 2)]
			lengths = [(p[1] - p[0]).getMag() for p in lines]
			#remove duplicates
			lengths = list(set(lengths))
			#remove the "diagonal" length
			lengths.remove(max(lengths))
			#if its a square (i.e. two identical lengths)
			if(len(lengths) == 1):
				return (lengths[0]**2)
			elif(len(lengths) == 2):
				return (lengths[0] * lengths[1])
			else:
				print "Something went wrong in Region.py area calculation"
				return 0
				
		elif(self.shape == 'circle'):
			return (np.pi*(self.args[0])**2)

		elif(self.shape == 'cylinder'):
			circumf = np.pi*2*self.args[0]
			l = self.args[1]
			return l*circumf

		return 0


		

	
	#return one of the two normal vectors
	#to the square plane region
	def getNormal(self):
		#plane defined by 3 points.
		if(self.shape == 'square'):
			p1 = self.args[1]
			p2 = self.args[2]
			p3 = self.args[3]
		elif(self.shape == 'circle'):
			p1 = Point.Point(self.args[0], 0, self.args[1], 0)
			p2 = Point.Point(0, 0, self.args[1], 0)
			p3 = Point.Point(self.args[0], 20, self.args[1], 0)
		elif(self.shape == 'cylinder'):
			return None

		#two vectors from 3 points
		v1 = p2 - p1
		v2 = p3 - p1

		#cross product gives a normal vector, then normalize 
		n1 = v2.cross(v1).normalize()
		n2 = v1.cross(v2).normalize()

		#in most cases one does not care which direction is
		#"outward" or "inward". For now, just return n1
		return n1

	#return the coefficients of the polynomial
	#Ax + By + Cz - D = 0
	#which define the plane. Returned in order
	#(A, B, C, D)
	def getPlaneCoefficients(self):
		n = self.getNormal()
		a = n.getX()
		b = n.getY()
		c = n.getZ()

		#shift plane from origin by one of the points
		#that defines it
		if(self.shape == 'square'):
			p = self.args[0]
		elif(self.shape == 'circle'):
			p = Point.Point(0, 0, self.args[1], 0)
		elif(self.shape == 'cylinder'):
			return None

		px = p.getX()
		py = p.getY()
		pz = p.getZ()

		d = a*px + b*py + c*pz
		return (a,b,c,d)

	def findIntersectionPoint(self, ray):
		if(len(ray.track) < 2):
			print "Ray was just a point! Please propagate before determining intersections"
			sys.exit()


		#computation for square or circle is different from cylinder
		if(self.shape == 'circle' or self.shape == 'square'):
			#plane equation coefficients
			(a, b, c, d) = self.getPlaneCoefficients()
			([x1, x2], [y1, y2], [z1, z2]) = ray.getCurrentLineParameters()
			
			#line is parametrized in t
			#point of intersection is found with the parameter
			epsilon = 0.000001
			if(abs((a*x2 + b*y2 + c*z2)) < epsilon):
				return None
			t = (d - a*x1 - b*y1 - c*z1)/(a*x2 + b*y2 + c*z2)
			#if t is negative, the intersection is in the past
			if(t < 0):
				return None
			x = x1 + x2*t
			y = y1 + y2*t
			z = z1 + z2*t
			return Point.Point(x,y,z,1)
		elif(self.shape == 'cylinder'):
			#for sure absorb the ray if it escaped
			#the boundary of the cylinder
			currentP = ray.getPos()
			if(currentP.getRho() > self.args[0] or currentP.getZ() < 0 or currentP.getZ() > self.args[1]):
				ray.setAbsorbed(True)


			([x1, x2], [y1, y2], [z1, z2]) = ray.getCurrentLineParameters()
			#find where it intersects the equation for the cylinder
			#x^2 + y^2 = r^2
			roots = np.roots([x2**2 + y2**2, 2*x1*x2 + 2*y1*y2, x1**2 + y1**2 - self.args[0]**2])
			#positive root is the intersection in the future
			t = max(roots)
			if(t < 0):
				return None
			x = x1 + x2*t
			y = y1 + y2*t
			z = z1 + z2*t
			return Point.Point(x, y, z, 1)

	#if the shape is a square,
	#return the bounds in phi
	def getPhiBounds(self):
		if(self.shape != 'square'):
			return None
		else:
			phi1 = self.args[0].getPhi()
			phi2 = None
			for s in self.args:
				if(phi1 - s.getPhi() != 0):
					phi2 = s.getPhi()

			phirange = [min([phi1, phi2]), max([phi1, phi2])]
			return phirange

	def getZBounds(self):
		if(self.shape != 'square'):
			return None
		else:
			zs = [f.getZ() for f in self.args]
			zdifs = [z[1] - z[0] for z in combinations(zs, 2)]
			return (min(zs) + max(zdifs))


	#checks if a given point is within the region boundaries
	def isPointInRegion(self, p):
		eps = 1e-6 #1 micron
		if(self.shape == 'circle'):
			rrange = [0, self.args[0]]
			zrange = [self.args[1] - eps, self.args[1] + eps]
			if(p.getRho() > rrange[0] and p.getRho() < rrange[1] \
				and p.getZ() > zrange[0] and p.getZ() < zrange[1]):
				return True
			else:
				return False

		elif(self.shape == 'cylinder'):
			rrange = [self.args[0] - eps, self.args[0] + eps]
			zrange = [0, self.args[1]]
			if(p.getRho() > rrange[0] and p.getRho() < rrange[1] \
				and p.getZ() > zrange[0] and p.getZ() < zrange[1]):
				return True
			else:
				return False

		elif(self.shape == 'square'):
			#check that it is within a phi range
			#and a z range. No constraint on Rho(r)

			#two pairs of points have z1 - z2 = 0
			#two pairs have z1 - z2 = K
			z1 = self.args[0].getZ()
			z2 = None
			for s in self.args:
				if(z1 - s.getZ() != 0):
					z2 = s.getZ()
			
			zrange = [min([z1,z1]), max([z1,z2])]
			#print "Point phi: " + str(p.getPhi())
			#print "ZRANGE = ",
			#print zrange

			#similar computation for phi range
			phi1 = self.args[0].getPhi()
			phi2 = None
			for s in self.args:
				if(phi1 - s.getPhi() != 0):
					phi2 = s.getPhi()

			phirange = [min([phi1, phi2]), max([phi1, phi2])]
			#print "PHIRANGE = ",
			#print phirange

			if(p.getZ() > zrange[0] and p.getZ() < zrange[1] and \
				p.getPhi() > phirange[0] and p.getPhi() < phirange[1]):
				return True
			else:
				return False



