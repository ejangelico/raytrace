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
		self.W = None
		self.H = None
		self.mirrorEnds = False
		self.type = None #describes the global geometry, "sphere", "cylinder", "rectangular"

	def getRegions(self):
		return self.regions

	def getL(self):
		return self.L

	def getR(self):
		return self.R




	#purely mirrored cylinder
	def mirrorCylinder(self, R, L):
		self.R = R
		self.L = L
		self.type = "cylinder"
		self.mirrorEnds = True

		cylReg = Region.Region('cylinder', 'mirror', [self.R, self.L])
		self.regions.append(cylReg)

		#do the end caps
		if(self.mirrorEnds):
			self.regions.append(Region.Region('circle', 'mirror', [self.R, -0.5*self.L]))
			self.regions.append(Region.Region('circle', 'mirror', [self.R, 0.5*self.L]))
		else:
			self.regions.append(Region.Region('circle', 'absorber', [self.R, -0.5*self.L]))
			self.regions.append(Region.Region('circle', 'absorber', [self.R, 0.5*self.L]))

	#purely mirrored sphere
	def mirrorSphere(self, R):
		self.R = R
		self.type = "sphere"
		self.regions.append(Region.Region('sphere', 'mirror', [self.R]))

	#purely mirrored rectangle
	def mirrorRect(self, l, w, h):
		self.W = w
		self.L = l
		self.H = h
		self.regions.append(Region.Region("rectVolume", 'mirror', [self.L, self.W, self.H]))




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
		incline = 0
		for a in z:
			for p in phi:
				self.makeCylLAPPD(p, a, incline, mat)

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
					self.makeCylLAPPD(p, z[i], incline, mats[0])
				else:
					self.makeCylLAPPD(p, z[i], incline, mats[1])


	#general detector with mirrors that are 
	#parallel with the LAPPDs

	#N = number of LAPPDs per column
	#n = number of columns
	def loadNonTilted(self, n, N, radius, mirrorEnds):
		#create a triangle pattern
		#two of the triangle edges are mirrors
		#the other is photocathode

		self.R = radius
		self.mirrorEnds = mirrorEnds


		l = 0.219964 #meters length of photocathode
		#adjust R based on minimum value
		#for the requested number "n"
		if(self.R <= l/(2.0*np.sin(np.pi/float(2.0*n)))):
			slop = 0.01 #1 cm of slop
			self.R = l/(2.0*np.sin(np.pi/float(2.0*n))) 

		#assume a cherenkov angle of 42
		effectiveR = 0.5*np.sqrt(4*self.R*self.R - l*l)
		sep = (2*effectiveR/np.tan(self.dtor(42))) - l
		self.L = N*(sep + l) + sep
		#cathode column mirror length are equal to sep
		#with a lower bound at 1/8 of lappd length
		if(sep < (1/8.0)*l):
			mirL = 0
		else:
			mirL = sep

		#calculate phi's
		phiInt = 360.0/n
		cathodePhis = [i*phiInt for i in range(n)]
		mirrorPhis = [(i*phiInt + phiInt/2.0) for i in range(n)]

		#calculate z's for cathode columns
		#based on there being N+1 separations
		#and N photocathodes
		zs = []
		for i in range(2*N + 1):
			if(i == 0):
				zs.append([0, 'mirror'])
			elif(i % 2 != 0):
				zs.append([zs[-1][0] + sep, 'cathode'])
			else:
				zs.append([zs[-1][0] + l, 'mirror'])

		#form photocathode columns
		for z in zs:
			for pcp in cathodePhis:
				if(z[1] == 'mirror'):
					if(mirL == 0):
						continue
					else:
						self.makeCylRectangle(pcp, z[0], 0, z[1], l, mirL)
				else:
					self.makeCylRectangle(pcp, z[0], 0, z[1], l, l)

		#form mirror columns
		for z in zs:
			for mp in mirrorPhis:
				if(z[1] == 'mirror'):
					if(mirL == 0):
						continue
					else:
						self.makeCylRectangle(mp, z[0], 0, 'mirror', l, mirL)
				else:
					self.makeCylRectangle(mp, z[0], 0, 'mirror', l, l)


		#make detector body
		cylReg = Region.Region('cylinder', 'absorber', [self.R, self.L])
		self.regions.append(cylReg)

		#do the end caps
		if(self.mirrorEnds):
			self.regions.append(Region.Region('circle', 'mirror', [self.R, 0]))
			self.regions.append(Region.Region('circle', 'mirror', [self.R, self.L]))
		else:
			self.regions.append(Region.Region('circle', 'absorber', [self.R, 0]))
			self.regions.append(Region.Region('circle', 'absorber', [self.R, self.L]))


		##Detector Header
		print "***Detector Header***"
		print str(n) + " columns of " + str(N) + " LAPPDs each "
		print "cylindrical volume length " + str(self.L)
		print "cylindrical volume radius " + str(self.R)
		print "water volume = " + str(np.pi*self.R*self.R*self.L*1000) + "L"
		print "separation of LAPPDs: " + str(sep)
		print "separating mirror length: " + str(mirL)
		f = self.coverageFractions()
		print "Total photocathode area = " + str(f[0]) + "m^2"
		print "PC Area/Mirror Area = " + str(f[2])
		print "PC Area/(Mirror + PC Area) = " + str(f[1])


	
	#a detector with tilted mirrors to reflect
	#photons right into the opposing LAPPD	
	#N = number of LAPPDs per column
	#n = number of columns
	#rsigma = 1sig uncertainty in R due to tilted mirros
	#rreduce = the reduction of the radius that reflection photons take
	#due to finite length of tilted mirrors
	def loadTilted(self, n, N, rsigma, rreduce, radius, mirrorEnds):


		self.R = radius
		self.mirrorEnds = mirrorEnds


		l = 0.219964 #meters length of photocathode
		#adjust R based on minimum value
		#for the requested number "n"
		if(self.R <= l/(2.0*np.sin(np.pi/float(2.0*n)))):
			slop = 0.01 #1 cm of slop
			self.R = l/(2.0*np.sin(np.pi/float(2.0*n))) 

		#assume a cherenkov angle of 42
		effectiveR = 0.5*np.sqrt(4*self.R*self.R - l*l)
		sep = (2*effectiveR/np.tan(self.dtor(42))) - l
		self.L = N*(sep + l) + sep
		#cathode column mirror length are equal to sep
		#with a lower bound at 1/8 of lappd length
		if(sep < (1/8.0)*l):
			mirL = 0
		else:
			mirL = sep

		#calculate phi's
		phiInt = 360.0/n
		cathodePhis = [i*phiInt for i in range(n)]
		mirrorPhis = [(i*phiInt + phiInt/2.0) for i in range(n)]

		#calculate z's for cathode columns
		#based on there being N+1 separations
		#and N photocathodes
		zs = []
		for i in range(2*N + 1):
			if(i == 0):
				zs.append([0, 'mirror'])
			elif(i % 2 != 0):
				zs.append([zs[-1][0] + sep, 'cathode'])
			else:
				zs.append([zs[-1][0] + l, 'mirror'])

		#form photocathode columns
		for z in zs:
			for pcp in cathodePhis:
				if(z[1] == 'mirror'):
					if(mirL == 0):
						continue
					else:
						self.makeCylRectangle(pcp, z[0], 0, z[1], l, mirL)
				else:
					self.makeCylRectangle(pcp, z[0], 0, z[1], l, l)


		#here is where the mirror folding
		#comes into play. One must now alter
		#the mirror attributes to be tilted and
		#short and have small separation, but only
		#for mirrors that oppose photocathode (i.e.,
		#mirrors with z[1] == 'cathode')
		
		#the given parameter rsigma is the "r" distance
		#from the center of exposed tilted mirror 
		sm = 2.0*rsigma/np.tan(self.dtor(24)) 	#separation of small mirrors based on 42 cherenkov angle
		mirrorWidth = rreduce/np.sin(self.dtor(24)) + rsigma/(np.sin(self.dtor(24)))
		tiltzs = np.linspace(0, l, int(l/sm))

		#form mirror columns
		for z in zs:
			for mp in mirrorPhis:
				if(z[1] == 'mirror'):
					if(mirL == 0):
						continue
					else:
						self.makeCylRectangle(mp, z[0], 0, 'mirror', l, mirL)
				else:
					#do a set of small tiled mirrors
					for tz in tiltzs:
						self.makeCylRectangle(mp, z[0] + tz, 24, 'mirror', l, mirrorWidth)


		#make detector body
		cylReg = Region.Region('cylinder', 'absorber', [self.R, self.L])
		self.regions.append(cylReg)

		#do the end caps
		if(self.mirrorEnds):
			self.regions.append(Region.Region('circle', 'mirror', [self.R, 0]))
			self.regions.append(Region.Region('circle', 'mirror', [self.R, self.L]))
		else:
			self.regions.append(Region.Region('circle', 'absorber', [self.R, 0]))
			self.regions.append(Region.Region('circle', 'absorber', [self.R, self.L]))


		##Detector Header
		print "***Detector Header***"
		print str(n) + " columns of " + str(N) + " LAPPDs each "
		print "cylindrical volume length " + str(self.L)
		print "cylindrical volume radius " + str(self.R)
		print "water volume = " + str(np.pi*self.R*self.R*self.L*1000) + "L"
		print "effective phi-slice radius " + str(effectiveR)
		print "separation of LAPPDs: " + str(sep)
		print "separating mirror length: " + str(mirL)
		print "separation of small mirrors: " + str(sm)
		print "number of small mirrors: " + str(len(tiltzs))
		print "effective bounce radius off tilted mirrors: " + str(effectiveR - rreduce)
		print "uncertainty on that bounce radius: " + str(rsigma)
		f = self.coverageFractions()
		print "Total photocathode area = " + str(f[0]) + "m^2"
		print "PC Area/Mirror Area = " + str(f[2])
		print "PC Area/(Mirror + PC Area) = " + str(f[1])


	#at this point I made a decision to make the length of the cylinder
	#an integer multiple of LAPPD lengths, with a mirror at both ends of 
	#each column

	#also, with V2, I've added tilted mirrors everywhere in the 
	#"send that photon back where it came from" design, for maximum
	#light collection

	#I changed direction before getting this to work
	def loadTiltedV2(self, n, N, rsigma, rreduce, radius, mirrorEnds):


		self.R = radius
		self.mirrorEnds = mirrorEnds


		l = 0.219964 #meters length of photocathode
		#adjust R based on minimum value
		#for the requested number "n"
		if(self.R <= l/(2.0*np.sin(np.pi/float(2.0*n)))):
			slop = 0.01 #1 cm of slop
			self.R = l/(2.0*np.sin(np.pi/float(2.0*n))) 

		#all mirror sections are the length of an LAPPD
		nsquares = 2*N + 1 	#number of "squares", either PC or Mirror, in each column
		self.L = l*nsquares

		#length of all mirrors is length of LAPPD
		mirL = l

		#calculate phi's
		phiInt = 360.0/n
		cathodePhis = [i*phiInt for i in range(n)]
		mirrorPhis = [(i*phiInt + phiInt/2.0) for i in range(n)]


		#calculate z's for cathode columns
		#based on there being N+1 separations
		#and N photocathodes
		zs = []
		for i in range(2*N + 1):
			if(i == 0):
				zs.append([0, 'mirror'])
			elif(i % 2 != 0):
				zs.append([zs[-1][0] + l, 'cathode'])
			else:
				zs.append([zs[-1][0] + l, 'mirror'])


		#calculate tilted mirror attributes for 
		#mirror sections opposite to PC columns
		#and mirror sections in PC columns

		#mirrors opposite photocathode
		sm = 2.0*rsigma/np.tan(self.dtor(24)) 	#separation of small mirrors based on 42 cherenkov angle
		mirrorWidth = rreduce/np.sin(self.dtor(24)) + rsigma/(np.sin(self.dtor(24)))
		tiltzs = np.linspace(0, l, int(l/sm))

		#mirrors opposite mirrors
		#in mirror columns
		sm2 = 2.0*rsigma/np.tan(self.dtor(90 - 42))
		mirrorWidth2 = rreduce/np.sin(self.dtor(90 - 42)) + rsigma/(np.sin(self.dtor(90 - 42)))
		tiltzs2 = np.linspace(0, l, int(l/sm2))



		#form photocathode columns
		for z in zs:
			for pcp in cathodePhis:
				if(z[1] == 'mirror'):
					if(mirL == 0):
						continue
					else:
						#do a set tilted mirrors for 
						#photocathode column angle
						for tz in tiltzs:
							self.makeCylRectangle(pcp, z[0] + tz, 24, 'mirror', l, mirrorWidth)
				else:
					self.makeCylRectangle(pcp, z[0], 0, z[1], l, l)

		#form mirror columns
		for z in zs:
			for mp in mirrorPhis:
				if(z[1] == 'mirror'):
					if(mirL == 0):
						continue
					else:
						#do a set of mirrors
						#for mirror column, high tilt
						for tz in tiltzs2:
							self.makeCylRectangle(mp, z[0] + tz, 90 - 42, 'mirror', l, mirrorWidth2)
				else:
					#do a set of small tiled mirrors
					for tz in tiltzs:
						self.makeCylRectangle(mp, z[0] + tz, 24, 'mirror', l, mirrorWidth)


		#make detector body
		cylReg = Region.Region('cylinder', 'absorber', [self.R, self.L])
		self.regions.append(cylReg)

		#do the end caps
		if(self.mirrorEnds):
			self.regions.append(Region.Region('circle', 'mirror', [self.R, 0]))
			self.regions.append(Region.Region('circle', 'mirror', [self.R, self.L]))
		else:
			self.regions.append(Region.Region('circle', 'absorber', [self.R, 0]))
			self.regions.append(Region.Region('circle', 'absorber', [self.R, self.L]))



		effectiveR = 0.5*np.sqrt(4*self.R*self.R - l*l)
		##Detector Header
		print "***Detector Header***"
		print str(n) + " columns of " + str(N) + " LAPPDs each "
		print "cylindrical volume length " + str(self.L)
		print "cylindrical volume radius " + str(self.R)
		print "water volume = " + str(np.pi*self.R*self.R*self.L*1000) + "L"
		print "effective phi-slice radius " + str(effectiveR)
		print "separation of low angle mirrors: " + str(sm)
		print "separation of high angle mirrors: " + str(sm2)
		print "number of low angle mirrors: " + str(len(tiltzs))
		print "number of high angle mirrors: " + str(len(tiltzs2))
		print "effective bounce radius off tilted mirrors: " + str(effectiveR - rreduce)
		print "uncertainty on that bounce radius: " + str(rsigma)
		f = self.coverageFractions()
		print "Total photocathode area = " + str(f[0]) + "m^2"
		print "PC Area/Mirror Area = " + str(f[2])
		print "PC Area/(Mirror + PC Area) = " + str(f[1])



	#makes an LAPPD sized photocathode
	#whose corner (closest to origin) is
	#at phi and z and r in cylindrical coordinates

	#the edge of the LAPPD is in contact
	#with the walls of the cylindrical detector
	#i.e. it makes a secant of the cylinder
	def makeCylLAPPD(self, phi, z, incline, mat):

		#basically need 4 points to define the 
		#square region. 
		lappd = 0.219964 #meters length of photocathode
		#find another edge that makes an intersection
		#with the cylinder wall, going in the +phi direction
		phicenter = phi
		delphi = 2*np.arcsin(lappd/(2*self.R))
		philow = phicenter - delphi/2.0
		phihi = phicenter + delphi/2.0

		#start with the square aligned along the x axis at z = 0
		p1 = Point.Point(lappd/2.0, 0, 0, 1)
		p2 = Point.Point(-lappd/2.0, 0, 0, 1)
		p3 = Point.Point(lappd/2.0, 0, lappd, 1)
		p4 = Point.Point(-lappd/2.0, 0, lappd, 1)

		#tilt it down by the incline only in the y direction
		ycor = Point.Point(0, lappd*np.sin(self.dtor(incline)), 0, 1)
		zcor = Point.Point(0, 0, lappd*(1 - np.cos(self.dtor(incline))), 1)
		p3 = p3 - zcor + ycor
		p4 = p4 - zcor + ycor

		#rotate all points about the z-axis, phi degrees
		rphi = self.dtor(phi)
		rotmat = np.matrix([[np.cos(rphi), -np.sin(rphi), 0], [np.sin(rphi), np.cos(rphi), 0], [0, 0, 1]])
		p1 = p1.rotate(rotmat)
		p2 = p2.rotate(rotmat)
		p3 = p3.rotate(rotmat)
		p4 = p4.rotate(rotmat)

		#shift the plane so that the lower z points
		#are touching the walls of the detector
		alpha = 2*np.arcsin(lappd/(4.0*self.R)) #secant angle 
		hype = self.R*np.sin(np.pi/2.0 - alpha)
		ycor = Point.Point(0, hype*np.sin(np.pi + self.dtor(phicenter + 90)), 0, 1)
		xcor = Point.Point(hype*np.cos(np.pi + self.dtor(phicenter + 90)), 0, 0, 1)
		#shift by z
		zcor = Point.Point(0, 0, z, 1)
		p1 = p1 + xcor + ycor + zcor
		p2 = p2 + xcor + ycor + zcor
		p3 = p3 + xcor + ycor + zcor
		p4 = p4 + xcor + ycor + zcor





		self.regions.append(Region.Region('square', mat, [p1, p2, p3, p4]))

	def makeCylRectangle(self, phi, z, incline, mat, length, width):

		#basically need 4 points to define the 
		#square region. 

		phicenter = phi
		delphi = 2*np.arcsin(length/(2*self.R))
		philow = phicenter - delphi/2.0
		phihi = phicenter + delphi/2.0

		#start with the square aligned along the x axis at z = 0
		p1 = Point.Point(length/2.0, 0, 0, 1)
		p2 = Point.Point(-length/2.0, 0, 0, 1)
		p3 = Point.Point(length/2.0, 0, width, 1)
		p4 = Point.Point(-length/2.0, 0, width, 1)

		#tilt it down by the incline only in the y direction
		ycor = Point.Point(0, width*np.sin(self.dtor(incline)), 0, 1)
		zcor = Point.Point(0, 0, width*(1 - np.cos(self.dtor(incline))), 1)
		p3 = p3 - zcor + ycor
		p4 = p4 - zcor + ycor

		#rotate all points about the z-axis, phi degrees
		rphi = self.dtor(phi)
		rotmat = np.matrix([[np.cos(rphi), -np.sin(rphi), 0], [np.sin(rphi), np.cos(rphi), 0], [0, 0, 1]])
		p1 = p1.rotate(rotmat)
		p2 = p2.rotate(rotmat)
		p3 = p3.rotate(rotmat)
		p4 = p4.rotate(rotmat)

		#shift the plane so that the lower z points
		#are touching the walls of the detector
		alpha = 2*np.arcsin(length/(4.0*self.R)) #secant angle 
		hype = self.R*np.sin(np.pi/2.0 - alpha)
		ycor = Point.Point(0, hype*np.sin(np.pi + self.dtor(phicenter + 90)), 0, 1)
		xcor = Point.Point(hype*np.cos(np.pi + self.dtor(phicenter + 90)), 0, 0, 1)
		#shift by z
		zcor = Point.Point(0, 0, z, 1)
		p1 = p1 + xcor + ycor + zcor
		p2 = p2 + xcor + ycor + zcor
		p3 = p3 + xcor + ycor + zcor
		p4 = p4 + xcor + ycor + zcor


		self.regions.append(Region.Region('square', mat, [p1, p2, p3, p4]))

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
		for reg in self.regions:
			sect = reg.findIntersectionPoint(ray)
			if(isinstance(sect, list)):
				for s in sect:
					if(s == None):
						continue
					elif(isinstance(s.x0, complex) or isinstance(s.x1, complex)):
						#point is outside cylindrical detector
						#the equation gives complex intersection point sometimes if outside of cylinder
						continue
					else:
						intersections.append([s, reg])
			else:
				if(sect == None):
					continue
				elif(isinstance(sect.x0, complex) or isinstance(sect.x1, complex)):
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


