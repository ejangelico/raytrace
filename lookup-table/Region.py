import Point

import sys
import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from itertools import combinations, product
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
		#sphere [radius]
		#rectVolume [x, y, z] dimensions
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
		elif(self.shape == 'rectVolume'):
			self.drawRectVolume(ax)
		elif(self.shape == 'sphere'):
			self.drawSphere(ax)


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

		nlines = 12.0
		length = self.args[1]
		dl = length/nlines
		dphi = 2*np.pi/nlines
		r = self.args[0]
		plotpoints = []
		angs = np.linspace(0, 2*np.pi, nlines, endpoint=True)
		for a in angs:
			plotpoints.append(Point.Point(r, a, -0.5*length, 0))
			plotpoints.append(Point.Point(r, a, 0.5*length, 0))
			x = [p.getX() for p in plotpoints]
			y = [p.getY() for p in plotpoints]
			z = [p.getZ() for p in plotpoints]
			ax.plot(x, y, z, color=self.getColor(), linestyle='-')
			plotpoints = []

		lens = np.linspace(-0.5*length, 0.5*length, nlines, endpoint=True)
		angs = np.linspace(0, 2*np.pi, 30.0, endpoint=True)
		for l in lens:
			for a in angs:
				plotpoints.append(Point.Point(r, a, l, 0))

			x = [p.getX() for p in plotpoints]
			y = [p.getY() for p in plotpoints]
			z = [p.getZ() for p in plotpoints]
			ax.plot(x, y, z, color=self.getColor())
			plotpoints = []


	def drawSphere(self, ax):
		u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
		x = self.args[0]*np.cos(u)*np.sin(v)
		y = self.args[0]*np.sin(u)*np.sin(v)
		z = self.args[0]*np.cos(v)
		ax.plot_wireframe(x, y, z, color=self.getColor())

	def drawRectVolume(self, ax):
		l = self.args[0]
		w = self.args[1]
		h = self.args[2]

		x = [-0.5*l, 0.5*l]
		y = [-0.5*w, 0.5*w]
		z = [-0.5*h, 0.5*h]
		xx, yy = np.meshgrid(x, y)
		for value in [0, 1]:
			ax.plot_wireframe(xx, yy, z[value], color=self.getColor())

		yy, zz = np.meshgrid(y, z)
		for value in [0, 1]:
			ax.plot_wireframe(x[value], yy, zz, color=self.getColor())

		xx, zz = np.meshgrid(x, z)
		for value in [0, 1]:
			ax.plot_wireframe(xx, y[value], zz, color=self.getColor())


	def getLengths(self):
		if(self.shape != 'square'):
			return None
		else:
			lines = [x for x in combinations(self.args, 2)]
			lengths = [(p[1] - p[0]).getMag() for p in lines]
			#remove duplicates
			return lengths

	def getArea(self):
		if(self.shape == 'square'):
			lines = [x for x in combinations(self.args, 2)]
			lengths = [(p[1] - p[0]).getMag() for p in lines]
			lengthCombs = [x for x in combinations(lengths, 2)]
			posAreas = [x[0]*x[1] for x in lengthCombs]
			posAreas = list(set(posAreas))
			#area is always the item that is next highest from the minimum
			#possible area
			posAreas = sorted(posAreas)
			#exactly square, then are is minimum
			if(len(posAreas) == 3):
				return posAreas[0]
			else:
				return posAreas[1]
			
				
		elif(self.shape == 'circle'):
			return (np.pi*(self.args[0])**2)

		elif(self.shape == 'cylinder'):
			circumf = np.pi*2*self.args[0]
			l = self.args[1]
			return l*circumf

		elif(self.shape == 'sphere'):
			return 4*np.pi*(self.args[0])**2

		elif(self.shape == 'rectVolume'):
			l = self.args[0]
			w = self.args[1]
			h = self.args[2]

			return 2*(w*l + h*l + h*w)

		return 0


		

	
	#return one of the two normal vectors
	#of the region at point "atP"
	def getNormal(self, atP):
		#plane defined by 3 points.
		if(self.shape == 'square'):
			p1 = self.args[1]
			p2 = self.args[2]
			p3 = self.args[3]
			#two vectors from 3 points
			v1 = p2 - p1
			v2 = p3 - p1

			#cross product gives a normal vector, then normalize 
			n1 = v2.cross(v1).normalize()

			#in most cases one does not care which direction is
			#"outward" or "inward". For now, just return n1
			return n1
		elif(self.shape == 'circle'):
			p1 = Point.Point(self.args[0], 0, self.args[1], 0)
			p2 = Point.Point(0, 0, self.args[1], 0)
			p3 = Point.Point(self.args[0], 20, self.args[1], 0)
			#two vectors from 3 points
			v1 = p2 - p1
			v2 = p3 - p1

			#cross product gives a normal vector, then normalize 
			n1 = v2.cross(v1).normalize()

			#in most cases one does not care which direction is
			#"outward" or "inward". For now, just return n1
			return n1
		elif(self.shape == 'cylinder'):
			n1 = Point.Point(-1*atP.getRho(), atP.getPhi(), 0, 0)
			return n1.convertToCoord(1).normalize()
		elif(self.shape == "sphere"):
			p1 = atP.convertToCoord(1)
			p2r = p1.getR() - 0.1*p1.getR()
			p2t = p1.getTh()
			p2p = p1.getPhi()
			cartP2 = Point.Point(p2r*np.cos(p2p)*np.sin(p2t), p2r*np.sin(p2p)*np.sin(p2t), p2r*np.cos(p2t), 1)
			n1 = (cartP2 - p1).normalize()
			return n1
		elif(self.shape == "rectVolume"):
			#make a bunch of square regions 
			#and find intersections with them
			l = self.args[0]
			w = self.args[1]
			h = self.args[2]

			x = [-0.5*l, 0.5*l]
			y = [-0.5*w, 0.5*w]
			z = [-0.5*h, 0.5*h]
			z0 = Point.Point(x[1], y[1], z[0], 1)
			z1 = Point.Point(x[0], y[1], z[0], 1)
			z2 = Point.Point(x[1], y[0], z[0], 1)
			z3 = Point.Point(x[0], y[0], z[0], 1)
			#z = zd plane
			zz0 = Point.Point(x[1], y[1], z[1], 1)
			zz1 = Point.Point(x[0], y[1], z[1], 1)
			zz2 = Point.Point(x[1], y[0], z[1], 1)
			zz3 = Point.Point(x[0], y[0], z[1], 1)


			#six faces define a rectangular prism
			faces = [[z0,z1,zz0,zz1], \
					[z0,z2,zz0,zz2], \
					[z2,z3,zz2,zz3], \
			 		[z1,z3,zz1,zz3], \
			 		[zz0,zz1,zz2,zz3], \
			 		[z0,z1,z2,z3]]


			faceRegs = [Region('square', self.mat, f) for f in faces]
			for f in faceRegs:
				if(f.isPointInRegion(atP)):
					n1 = f.getNormal(atP)
					return n1

			return None

			


			






		

	#return the coefficients of the polynomial
	#Ax + By + Cz - D = 0
	#which define the plane. Returned in order
	#(A, B, C, D)
	def getPlaneCoefficients(self):
		
		#shift plane from origin by one of the points
		#that defines it
		if(self.shape == 'square'):
			p = self.args[0]
			n = self.getNormal(p)
		elif(self.shape == 'circle'):
			p = Point.Point(0, 0, self.args[1], 0)
			n = self.getNormal(p)
		else:
			return None


		a = n.getX()
		b = n.getY()
		c = n.getZ()

		px = p.getX()
		py = p.getY()
		pz = p.getZ()

		d = a*px + b*py + c*pz
		return (a,b,c,d)

	def findIntersectionPoint(self, ray):

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
			if(currentP.getRho() > (self.args[0] + 1e-8) or currentP.getZ() < -0.5*self.args[1] or currentP.getZ() > 0.5*self.args[1]):
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
			sectpoint = Point.Point(x, y, z, 1)
			return Point.Point(x, y, z, 1)

		elif(self.shape == 'sphere'):
			([x1, x2], [y1, y2], [z1, z2]) = ray.getCurrentLineParameters()
			#find where it intersects the equation 
			#x^2 + y^2 + z^2 = r^2
			roots = np.roots([x2**2 + y2**2 + z2**2, 2*x1*x2 + 2*y1*y2 + 2*z1*z2, x1**2 + y1**2 + z1**2 - self.args[0]**2])
			t = max(roots)
			if(t < 0):
				return None

			x = x1 + x2*t
			y = y1 + y2*t
			z = z1 + z2*t
			return Point.Point(x, y, z, 1)

		elif(self.shape == 'rectVolume'):
			#make a bunch of square regions 
			#and find intersections with them
			l = self.args[0]
			w = self.args[1]
			h = self.args[2]

			x = [-0.5*l, 0.5*l]
			y = [-0.5*w, 0.5*w]
			z = [-0.5*h, 0.5*h]
			z0 = Point.Point(x[1], y[1], z[0], 1)
			z1 = Point.Point(x[0], y[1], z[0], 1)
			z2 = Point.Point(x[1], y[0], z[0], 1)
			z3 = Point.Point(x[0], y[0], z[0], 1)
			#z = zd plane
			zz0 = Point.Point(x[1], y[1], z[1], 1)
			zz1 = Point.Point(x[0], y[1], z[1], 1)
			zz2 = Point.Point(x[1], y[0], z[1], 1)
			zz3 = Point.Point(x[0], y[0], z[1], 1)


			#six faces define a rectangular prism
			faces = [[z0,z1,zz0,zz1], \
					[z0,z2,zz0,zz2], \
					[z2,z3,zz2,zz3], \
			 		[z1,z3,zz1,zz3], \
			 		[zz0,zz1,zz2,zz3], \
			 		[z0,z1,z2,z3]]


			faceRegs = [Region('square', self.mat, f) for f in faces]

			intpoints = [frg.findIntersectionPoint(ray) for frg in faceRegs]
			return intpoints


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
			zrange = [-0.5*self.args[1], 0.5*self.args[1]]
			if(p.getRho() > rrange[0] and p.getRho() < rrange[1] \
				and p.getZ() > zrange[0] and p.getZ() < zrange[1]):
				return True
			else:
				return False

		elif(self.shape == 'square'):
			ps = [self.args[i].convertToCoord(1) for i in range(len(self.args) - 1)]
			#check if three of these points, plus the point in question
			#lies in the plane
			ps.append(p)
			M = [[_.getX() for _ in ps], \
				[_.getY() for _ in ps], \
				[_.getZ() for _ in ps], \
				[1,1,1,1]]

			#if determinant of this matrix is zero,
			#then the four points are coplanar
			if(abs(np.linalg.det(M)) > 1e-8):
				return False

			#now check if it is within the square
			#boundaries by getting max and min on
			#each x, y, and z for all points
			ps = [self.args[i].convertToCoord(1) for i in range(len(self.args))]
			carts = [[_.getX() for _ in ps], [_.getY() for _ in ps], [_.getZ() for _ in ps]]
			ranges = [[max(_) + 1e-8, min(_) - 1e-8] for _ in carts]
			if(ranges[0][0] >= p.getX() >= ranges[0][1] and \
				ranges[1][0] >= p.getY() >= ranges[1][1] and \
				ranges[2][0] >= p.getZ() >= ranges[2][1]):
				return True
			else:
				print "THIS is happening, check it out while you have the geometry"
				return False






		elif(self.shape == "sphere"):
			return True

		elif(self.shape == "rectVolume"):
			return True


