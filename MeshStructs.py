
import numpy as np

DIMENSION = 2

nNodes = -1
nElem = -1

ElemVertInds = None
VertCoords = None
groupNames = None
groupMembers = None

Areas = None
c0s = None
c1s = None
c2s = None

basisCoeffs = None #[element, xyc, vertex]

bodyNodeNormals = None

funcVolumes = None

""" Mesh Loading Routines """


def parseMesh(filename):


	global nNodes
	global nElem
	global ElemVertInds
	global VertCoords
	global Areas
	global c0s
	global c1s
	global c2s
	global basisCoeffs
	global groupNames
	global groupMembers
	global funcVolumes

	meshFile = open(filename,'r')
	
	#chop first line
	line = meshFile.readline()

	#get how many elements there are
	line = meshFile.readline().split()
	nElem = int(line[-1])

	#get indices of element vertices
	ElemVertInds = np.zeros([nElem,DIMENSION+1], dtype=int)
	for i in range(nElem):
		line = meshFile.readline().split()
		for j in range(DIMENSION+1):
			ElemVertInds[i,j] = int(line[j+1])


	#get how many vertices there are
	line = meshFile.readline().split()
	nNodes = int(line[-1])

	#get indices of cell vertices
	VertCoords = np.zeros([nNodes,DIMENSION])
	for i in range(nNodes):
		line = meshFile.readline().split()
		for j in range(DIMENSION):
			VertCoords[i,j] = float(line[j])


	# remove unconnected nodes
	fromInd = []
	for i in range(nNodes):
		if i in ElemVertInds:
			fromInd.append(i)
			ElemVertInds[ElemVertInds == i] = len(fromInd)-1
	VertCoords = VertCoords[fromInd,:]
	nNodes = VertCoords.shape[0]
	
	# PRINT VTK MESHFILE
	makeGridFile(ElemVertInds,VertCoords,"meshFile.vtk")


	#get how many groups there are
	line = meshFile.readline().split()
	nGroups = int(line[-1])


	groupSizes = np.zeros(nGroups, dtype=int)
	groupNames = ["" for _ in range(nGroups)]
	groupMembers = [[] for _ in range(nGroups)]

	for i in range(nGroups):
		#get group name
		line = meshFile.readline().split()
		groupNames[i] = line[-1]
		#get group size
		line = meshFile.readline().split()
		groupSizes[i] = int(line[-1])
		#get group data
		for j in range(groupSizes[i]):
			line = meshFile.readline().split()
			for k in range(2):
				ind = int(line[-(1+k)])
				if ind in fromInd: # must adjust for deleted nodes
					ind = fromInd.index(ind)
					
					if ind not in groupMembers[i]:
						groupMembers[i].append(ind)
		

	meshFile.close()
	
	# Computing element properties
	c0s = VertCoords[ ElemVertInds[:,0], : ]
	c1s = VertCoords[ ElemVertInds[:,1], : ]
	c2s = VertCoords[ ElemVertInds[:,2], : ]
	
	Areas = getElemAreas(c0s, c1s, c2s)
	
	funcVolumes = basisFuncVolumes()
	
	basisCoeffs = get2DBasisCoeffs(c0s, c1s, c2s) #[element, xyc, vertex]
	
	getBodyOutwardNormals()
# END PARSEMESH


def makeGridFile(vertexIndices, vertexCoords, filename):
	numVertices = vertexCoords.shape[0]
	numCells = vertexIndices.shape[0]
	
	vertices3D = np.zeros([numVertices,3])
	vertices3D[:,:2] = vertexCoords[:,:]
	
	cellData = 3*np.ones([numCells,4], dtype=int)
	cellData[:,1:] = vertexIndices[:,:]
	
	f = open(filename,'w')
	
	#header
	f.write("# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\n")
	
	#print coordinates of vertices
	f.write("POINTS " + str(numVertices) + " float\n")
	for i in range(numVertices):
		f.write(str(vertices3D[i,0]) + " " + str(vertices3D[i,1]) + " " + str(vertices3D[i,2]) + "\n" )

	
	#print coordinates of vertices
	f.write("CELLS " + str(numCells) + " " + str(4*numCells) + "\n")
	for i in range(numCells):
		f.write(str(cellData[i,0]) + " " + str(cellData[i,1]) + " " + str(cellData[i,2]) \
				+ " " + str(cellData[i,3]) + "\n" )
	
	#cell types
	f.write("CELL_TYPES " + str(numCells) + "\n")
	for i in range(numCells):
		f.write(str(5)+"\n")
	
	f.close()
#END OF makeGridFile


def printResults(rho,v1,v2,en, tag):
	from FluxFunctions import pressure
	
	P = pressure(rho,en)
	
	
	filename = "results.vtk." + str(tag)
	f = open(filename,'w')
	
	
	f.write( "POINT_DATA " + str(nNodes) + "\n" )
	
	#Scalar values
	f.write("SCALARS rho float\nLOOKUP_TABLE default\n")
	for i in range(nNodes):
		f.write( str(rho[i]) + "\n" )


	f.write("SCALARS e float\nLOOKUP_TABLE default\n")
	for i in range(nNodes):
		f.write(str(en[i]) + "\n")


	f.write("SCALARS T float\nLOOKUP_TABLE default\n")
	for i in range(nNodes):
		f.write(str(en[i]/718.0) + "\n")

	
	f.write("SCALARS P float\nLOOKUP_TABLE default\n")
	for i in range(nNodes):
		f.write(str(P[i]) + "\n")

	
	#Vector values
	f.write("VECTORS V float\n")
	for i in range(nNodes):
	
		f.write(str(v1[i]) + " " + str(v2[i]) + " " + str(0.0) + "\n")
	
	f.close()


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
	allMats[:,:,0] = np.column_stack((c0s[:,0],c1s[:,0],c2s[:,0]))
	allMats[:,:,1] = np.column_stack((c0s[:,1],c1s[:,1],c2s[:,1]))
	allMats[:,:,2] = np.ones([nElem,3])

	rhs = np.zeros([nElem,3,3])
	for i in range(3):
		rhs[:,i,i] = 1.0
	return np.linalg.solve(allMats, rhs) #[element, xyc, vertex]

""" End of Mesh Loading Routines """



""" 2D integration Routines """

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
		


def getIntegPtsWt(verts,levels):
	subTris = subdivideTriangle(verts,levels)
	
	evalPts = []
	wght = triArea(subTris[0]) # all triangles will have the same weight
	
	for tri in subTris:
		evalPts.append( np.mean(tri, axis=0) )
	
	return np.array(evalPts), wght


def integrateHatFunction(elem, vert, nlvls): # 2 levels should do it
	# vert is 0, 1, or 2
	coords = np.stack([c0s[elem,:],c1s[elem,:],c2s[elem,:]])
	intPts, w = getIntegPtsWt( coords, nlvls )
	
	return w*np.sum( intPts.dot(basisCoeffs[elem,:DIMENSION,vert]) + basisCoeffs[elem,DIMENSION,vert] )


def basisFuncVolumes():
	volumes = np.zeros( nNodes )
	
	for elem in range(nElem):
		for i in range(DIMENSION+1):
			volumes[ ElemVertInds[elem,i] ] += Areas[elem]/3.0
	
	return volumes


""" End 2D Integration Routines """

""" Sparse Matrix Construction Routines and Structures """

from scipy.sparse import coo_matrix, csr_matrix
class spMatBuilder:
	def __init__(self, N):
		self.N = N
		self.dat = []
		self.row = []
		self.col = []
		
		self.coo = None
		
	def addEntry(self,i,j,val):
		
		self.row.append(i)
		self.col.append(j)
		self.dat.append(val)
		
	def makeCOO(self):
			self.coo = coo_matrix( (self.dat, (self.row,self.col)), shape=(self.N,self.N) )
			self.coo.sum_duplicates()
			
			self.dat = self.coo.data.tolist()
			self.row = self.coo.row.tolist()
			self.col = self.coo.col.tolist()
	
	def getDense(self):
		if self.coo is None:
			self.makeCOO()
		return self.coo.toarray()
	
	def getSparse(self):
		if self.coo is None:
			self.makeCOO()
		return self.coo.tocsr()
	
	def getFullSparse(self):
		if self.coo is None:
			self.makeCOO()
		
		# 4 is b/c 4 variables
		
		allDats = [self.coo.data]
		allRows = [self.coo.row]
		allCols = [self.coo.col]
		
		for i in range(1,4):
			allDats.append( allDats[-1] )
			allRows.append( allRows[-1] + self.N )
			allCols.append( allCols[-1] + self.N )
		
		allDats = np.concatenate( allDats )
		allRows = np.concatenate( allRows )
		allCols = np.concatenate( allCols )
		
		return coo_matrix( (allDats, (allRows,allCols)), shape=(4*self.N,4*self.N) ).tocsr()


# VERIFIED UP TO HERE

def makeStiffnessMatrix():
	
	builder = spMatBuilder(nNodes)
	
	for elem in range(nElem):
		for i in range(DIMENSION+1):
			for j in range(DIMENSION+1):
				vert1 = ElemVertInds[elem,i]
				vert2 = ElemVertInds[elem,j]
				# TODO consider normalizing gradients to emphasize averaging term, not true Laplacian
				grad1 = basisCoeffs[elem,:DIMENSION,i] 
				grad2 = basisCoeffs[elem,:DIMENSION,j]
				
				val = Areas[elem]*grad1.dot(grad2)
				builder.addEntry(vert1,vert2,val)

	
	return builder


def makeMassMatrix():
	
	builder = spMatBuilder(nNodes)
	
	for elem in range(nElem):
		for i in range(DIMENSION+1):
			for j in range(DIMENSION+1):
				vert1 = ElemVertInds[elem,i]
				vert2 = ElemVertInds[elem,j]
				
				val = 1.0
				if vert1 == vert2:
					val = 2.0
				val = Areas[elem]*val/12.0
				builder.addEntry(vert1,vert2,val)

	
	return builder


def makeAvgMatrix():
	builder = makeMassMatrix()
	builder.makeCOO()
	
	for i in range(len(builder.dat)):
		builder.dat[i] /= funcVolumes[builder.row[i]]
	
	builder.makeCOO()
	
	return builder


def makeMixedMatrices():
	
	builder1 = spMatBuilder(nNodes)
	builder2 = spMatBuilder(nNodes)
	
	for elem in range(nElem):
		for i in range(DIMENSION+1):
			for j in range(DIMENSION+1):
				vert1 = ElemVertInds[elem,i]
				vert2 = ElemVertInds[elem,j]
				
				grad1 = basisCoeffs[elem,:DIMENSION,i] 
				
				val1 = Areas[elem]*grad1[0]/3.0
				val2 = Areas[elem]*grad1[1]/3.0
				
				builder1.addEntry(vert1,vert2,val1)
				builder2.addEntry(vert1,vert2,val2)

	# Boundary terms at the surface of the object
	GID = groupNames.index( 'Body' )	
	normal = np.zeros([DIMENSION])
	for elem in range(nElem):
		# see which nodes belong to the body
		bodyNodes = [ ElemVertInds[elem,i] in groupMembers[GID] for i in range(DIMENSION+1) ]
		
		# if a face of the element is a piece of the body
		if np.sum(bodyNodes) == DIMENSION:
			notBodyNodes = [not val for val in bodyNodes]
		
			vert1 = ElemVertInds[elem,bodyNodes][0]
			vert2 = ElemVertInds[elem,bodyNodes][1]
			nodeCoords = VertCoords[ ElemVertInds[elem,bodyNodes], : ]
			
			tangent = nodeCoords[1,:] - nodeCoords[0,:]
			dist = np.linalg.norm(tangent)
			
			outward = VertCoords[ ElemVertInds[elem,notBodyNodes], : ] - nodeCoords[0,:]
			
			# Compute normal direction and scale to unit vector
			normal[0] = tangent[1]
			normal[1] = -tangent[0]
			# Make sure is inward and unit
			normal *= -np.sign( outward.dot(normal) )
			normal /= np.linalg.norm(normal)
			
			builder1.addEntry(vert1,vert1, -(normal[0]*dist)/3.0 )
			builder1.addEntry(vert2,vert2, -(normal[0]*dist)/3.0 )
			builder1.addEntry(vert1,vert2, -(normal[0]*dist)/6.0 )
			builder1.addEntry(vert2,vert1, -(normal[0]*dist)/6.0 )
			
			builder2.addEntry(vert1,vert1, -(normal[1]*dist)/3.0 )
			builder2.addEntry(vert2,vert2, -(normal[1]*dist)/3.0 )
			builder2.addEntry(vert1,vert2, -(normal[1]*dist)/6.0 )
			builder2.addEntry(vert2,vert1, -(normal[1]*dist)/6.0 )

	
	return builder1, builder2

""" End Sparse Matrix Construction Routines and Structures """

def bodyConnectivity():
	GID = groupNames.index( 'Body' )

	neighbors = -np.ones( [len(groupMembers[GID]), DIMENSION] )
	
	for j in range(len(groupMembers[GID])):
		mmbr = groupMembers[GID][j]
		nhbr = 0
		elems = []
		for i in range(DIMENSION+1):
			elems.extend( np.where(ElemVertInds[:,i] == mmbr) )
		
		for elem in elems:
			for i in range(DIMENSION+1):
				if (ElemVertInds[elem,i] in groupMembers[GID]) and (ElemVertInds[elem,i] != mmbr) and (ElemVertInds[elem,i] != neighbors[j,0]) :
					neighbors[j,nhbr] = ElemVertInds[elem,i]
					nhbr += 1
			if nhbr >= 2:
				break
	return neighbors


def getBodyOutwardNormals():
	global bodyNodeNormals

	GID = groupNames.index( 'Body' )
	
	bodyNodeNormals = np.zeros([len(groupMembers[GID]),DIMENSION])
	
	normal = np.zeros([DIMENSION])
	
	for elem in range(nElem):
		# see which nodes belong to the body
		bodyNodes = [ ElemVertInds[elem,i] in groupMembers[GID] for i in range(DIMENSION+1) ]
		
		# if a face of the element is a piece of the body
		if np.sum(bodyNodes) == DIMENSION:
			notBodyNodes = [not val for val in bodyNodes]
		
			nodeCoords = VertCoords[ ElemVertInds[elem,bodyNodes], : ]
			
			tangent = nodeCoords[1,:] - nodeCoords[0,:]
			outward = VertCoords[ ElemVertInds[elem,notBodyNodes], : ] - nodeCoords[0,:]
			
			# Compute normal direction and scale to unit vector
			normal[0] = tangent[1]
			normal[1] = -tangent[0]
			# Make sure is outward
			normal *= np.sign( outward.dot(normal) )

			
			globalInds = ElemVertInds[elem,bodyNodes]
			for ind in globalInds:
				groupInd = groupMembers[GID].index(ind)
				bodyNodeNormals[groupInd,:] += normal
	
	# Take the average and normalize		
	bodyNodeNormals /= 2.0
	bodyNodeNormals /= np.sqrt(np.sum(bodyNodeNormals**2,1)).reshape([-1,1])
































