import numpy as np
import sys

#Vectors and Points can be same data structure
class Point:
	#pass three coordinates
	#point class is structured to
	#support cartesian or cylindrical coords
	def __init__(self, x0, x1, z, cart):
		self.x0 = x0 #x, rho, R
		self.x1 = x1 #y, phi, th
		self.z = z 	 #z, z, phi

		#is this point in cartesian or cylindrical? (True = Cartesian)
		self.cart = cart 

	def __mul__(self, b):
		return (self.getX() * b.getX()) + (self.getY() * b.getY()) + (self.getZ() * b.getZ())

	#scalar multiplication
	def scale(self, b):
		return Point(b*self.getX(), b*self.getY(), b*self.getZ(), 1)

	def __add__(self, b):
		#just return in cartesian

		return Point(self.getX() + b.getX(),self.getY() + b.getY(),self.getZ() + b.getZ(),1)  

	def __sub__(self, b):
		return (self + b.scale(-1))


	def __str__(self):
		return "(" + str(self.getX()) + ", " + str(self.getY()) + ", " + str(self.getZ()) + ")"

	def getMag(self):
		return np.sqrt(self.getX()**2 + self.getY()**2 + self.getZ()**2)

	def cross(self, b):
		#convert both vectors to cartesian
		x = [self.getX(), self.getY(), self.getZ()]
		y = [b.getX(), b.getY(), b.getZ()]
		xy = np.cross(x, y)
		return Point(xy[0], xy[1], xy[2], 1)

	def getCart(self):
		return self.cart

	def getX(self):
		a = self.convertToCoord(1)
		return a.x0

	def getY(self):
		a = self.convertToCoord(1)
		return a.x1


	#cylindrical coords
	def getPhi(self):
		a = self.convertToCoord(0)
		return a.x1

	def getRho(self):
		a = self.convertToCoord(0)
		return a.x0

	def getZ(self):
		return self.z

	#spherical coords
	def getR(self):
		a = self.convertToCoord(1)
		return (np.sqrt(a.x0**2 + a.x1**2 + a.z**2))

	def getTh(self):
		a = self.convertToCoord(1)
		r = a.getR()
		z = a.z
		if(z == 0):
			return np.pi/2.0
		elif(z/r >1 or z/r < -1):
			print "ERROR in calculation spherical theta"
			return 0
		else:
			return np.arccos(z/r)



	#checks if this point is effectively
	#the zero point or zero vector
	#accounting for floating point errors
	def isFloatZero(self):
		#check within 1 nm
		nm = 1e-9
		if(abs(self.getX()) < nm and \
		abs(self.getY()) < nm and \
		abs(self.getZ()) < nm ):
			return True
		else:
			return False


	def rotateIntoReferenceVector(self, refVec):
		R = self.getRotationMatrix(refVec)
		return self.rotate(R)

	#takes a vector "self" and returns a rotation martix
	#so that it points in the "newdirect"ion
	def getRotationMatrix(self, newdirect):
		normedSelf = self.normalize()
		normedDirect = newdirect.normalize()
		v = normedSelf.cross(normedDirect)
		s = v.getMag()
		c = normedSelf*normedDirect

		#happens if vectors are back to back, or the same
		epsilon = 0.0000001
		if(abs(c + 1) < epsilon):
			return np.matrix([[1,0,0],[0,1,0],[0,0,1]])

		ident = np.matrix([[1,0,0],[0,1,0],[0,0,1]])
		#skew symmetric cross-product matrix of v
		skew = np.matrix([[0, -1*v.getZ(), v.getY()],[v.getZ(), 0, -1*v.getX()],[-1*v.getY(), v.getX(), 0]]) 
		skew2 = skew*skew
		#is undefined if the two vectors are back to back
		rotationFactor = (1/(1 + c))
		R = ident + skew + skew2*rotationFactor
		return R


	#rotates the vector given a matrix
	#basically just matrix multiplication 
	#like Av = v'
	def rotate(self, rotMat):
		#only do rotations in cartesian

		#turn into a "numpy matrix column vector"
		v = np.matrix([[self.getX()],[self.getY()],[self.getZ()]])

		#expects rotMat to be a numpy matrix
		if(type(rotMat) != type(v)):
			print "Please make rotMat a np.matrix!"
			print "Tried: " + str(rotMat)
			print "STOPPING PROGRAM"
			sys.exit()

		newv = rotMat*v
		return Point(newv.item(0), newv.item(1), newv.item(2), 1)

	def normalize(self):
		mag = self.getMag()
		if(mag == 0):
			print "Tried to normalize the zero vector!"
			sys.exit()
		if(self.cart == 1):
			return self.scale(1.0/mag)
		else:
			#normalize in cartesian
			mag = self.getMag()
			b = self.convertToCoord(1)
			bNormed = b.scale(1.0/mag)
			#return in cylindrical
			return bNormed.convertToCoord(0)



	#convert to "coord"
	#1: Cartesian
	#0: cylindrical
	def convertToCoord(self, coord):
		#if it is the current coord system
		if(coord == self.cart):
			return self 
		elif(coord == 0):
			rho = np.sqrt(self.x0**2 + self.x1**2)
			phi = np.arctan2(self.x1, self.x0)
			if(phi < 0):
				phi = phi + 2*np.pi
			z = self.z
			return Point(rho, phi, z, coord)

		elif(coord == 1):
			x = self.x0 * np.cos(self.x1)
			y = self.x0 * np.sin(self.x1)
			z = self.z
			return Point(x, y, z, coord)


		else:
			print "Invalid input for coordinate frame, use 0 or 1: returning self"
			return self


