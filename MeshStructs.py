
import numpy as np


# ElemVertInds
# VertCoords
# if DIM == 2
#c0s = VertCoords[ ElemVertInds[:,0], : ]
#c1s = VertCoords[ ElemVertInds[:,1], : ]
#c2s = VertCoords[ ElemVertInds[:,2], : ]

#Areas = getElemAreas(c0s, c1s, c2s)

#basisCoeffs = get2DBasisCoeffs(c0s, c1s, c2s) #[element, xyc, vertex]

# form stiffness matrix

# form mass matrix

# form mixed matrix

def parseMesh(filename):
#filename = "FullCircle.su2"
	meshFile = open(filename,'r')
	
	#chop first line
	line = meshFile.readline()
	
	#get how many elements there are
	line = meshFile.readline().split()
	nElem = int(line[-1])
	
	#get indices of element vertices
	ElemVertInds = np.zeros([nElem,3]).astype(int)
	for i in range(nElem):
		line = meshFile.readline().split()
		for j in range(3):
			ElemVertInds[i,j] = int(line[j+1])
	
	
	#get how many vertices there are
	line = meshFile.readline().split()
	nNodes = int(line[-1])
	
	#get indices of cell vertices
	VertCoords = np.zeros([nNodes,2])
	for i in range(nNodes):
		line = meshFile.readline().split()
		for j in range(2):
			cellVertexCoords[i,j] = float(line[j])

	# TODO
#	makeGridFile(ElemVertInds,VertCoords,"meshFile.vtk")

	#get how many groups there are
	line = meshFile.readline().split()
	nGroups = int(line[-1])
	
	
	groupSizes = np.zeros(nGroups).astype(int)
	groupNames = ["" for _ in range(nGroups)]
	groupMembers = [[] for _ in range(nGroups)]
	
	for i in range(numGroups):
		#get group name
		line = meshFile.readline().split()
		groupNames[i] = line[-1]
		#get group size
		line = meshFile.readline().split()
		groupSizes[i] = int(line[-1])
		#get group data - not storing
		for j in range(groupSizes[i]):
			line = meshFile.readline().split()
			for k in range(2):
				ind = int(line[-(1+k)])
				if ind not in groupMembers[i]:
					groupMembers[i].append(ind)
			

	meshFile.close()




def triArea(verts):
	# verts is 3 by 2
	
	v1 = verts[0,:] - verts[1,:]
	v2 = verts[2,:] - verts[1,:]
	
	return np.abs( v1[0]*v2[1] - v1[1]*v2[0] )/2.0


def getElemAreas(c0s, c1s, c2s):
	V1s = c0s - c1s
	V2s = c2s - c1s

	return np.abs( V1s[:,0]*V2s[:,1] - V1s[:,1]*V2s[:,0] )/2.0


def get2DBasisCoeffs(c0s, c1s, c2s):
	# Computing the coefficients for the basis functions
	allMats = np.zeros([nElem,3,3])
	allMats[:,:,0] = np.concatenate((c0s[:,0],c1s[:,0],c2s[:,0]),1)
	allMats[:,:,1] = np.concatenate((c0s[:,1],c1s[:,1],c2s[:,1]),1)
	allMats[:,:,2] = np.ones([nElem,3])

	rhs = np.zeros([nElem,3,3])
	for i in range(nElem):
		rhs[i,:,:] = np.eye([3,3])
	return np.linalg.solve(allMats, rhs) #[element, xyc, vertex]



def subdivideTriangle(verts,levels):
	# verts is 3 by 2
	
	prevList = [verts]
	
	for _ in range(levels):
		newList = []
		for tri in prevList:
			mids = np.zeros([3,2])
			for i in range(3):
				pt1 = np.mod(i+1,3)
				pt2 = np.mod(i+2,3)
				mids[i,:] = (tri[pt1,:] + tri[pt2,:])/2.0

			newList.append( mids.copy() )

			for i in range(3):
				pt1 = np.mod(i+1,3)
				pt2 = np.mod(i+2,3)

				newTri = np.zeros([3,2])
				newTri[0,:] = tri[i,:]
				newTri[1,:] = mids[pt1,:]
				newTri[2,:] = mids[pt2,:]
				newList.append( newTri.copy() )

		prevList = newList.copy()

	return prevList
		

def getIntegPtsWts(verts,levels):
	subTris = subdivideTriangle(verts,levels)
	
	evalPts = []
	wghts = []
	
	for tri in subTris:
		evalPts.append( np.mean(tri, axis=0) )
		wghts.append( triArea(tri) )

	return evalPts, wghts



def keyedList:
	def __init__(self, ind=-1, content=[]):
		self.ind = ind
		self.content = content


def keyedValue:
	def __init__(self, ind=-1, val=0.0):
		self.ind = ind
		self.val = val

# list.sort(key=getListKey)
def getListKey(theList):
	return theList.ind

# list.sort(key=getValueKeyA)
def getValueKey(theValue):
	return theValue.ind


def spMatBuilder:
	def __init__(self):
		
	def addEntry(self,i,j,val):
		
		done = False
		for row in range(len(self.dat)):
			if self.dat[row].ind == i:
				# See if the column index exits
				for col in range(len(self.dat[row].content)):
					# if it exist, add the value to the entry
					if self.dat[row].content[col].ind == j:
						self.dat[row].content[col].val += val
						done = True
						break
				# if it doesn't exist, insert it
				if not done:
					self.dat[row].content.append( keyedValue(j,val) )
					done = True
					break
		if not done:
			self.dat.append( keyedList(i, [keyedValue(j,val)] ) )

