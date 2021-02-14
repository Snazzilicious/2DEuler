
import numpy as np
from scipy.sparse.linalg import spsolve
import sys
import importlib

if 'MeshStructs' in sys.modules:
	import MeshStructs as MS
	importlib.reload(MS)
else:
	import MeshStructs as MS

if 'FluxFunctions' in sys.modules:
	import FluxFunctions as FF
	importlib.reload(FF)
else:
	import FluxFunctions as FF


""" Set Control Variables """
visc = 1.0
numIts = 10
numSubIts = 10
NUM_VARS = 4

# Pick a mesh file
filename = "Mesh/2DCircle.su2"

# Load in the Mesh
print("Loading the Mesh", flush=True)
MS.parseMesh(filename)


""" Set free stream conditions """
M_inf = 1.25
AoA = 0*np.pi/180.0

# don't change these
Vinf = np.array([np.cos(AoA), np.sin(AoA)])
en_inf = ( 1.0 / ( FF.gamma*(FF.gamma-1)*M_inf*M_inf ) )



# Shortucts to Variables
rhoIndices = np.arange(MS.nNodes).astype(int)
v1Indices=np.arange(MS.nNodes).astype(int) + MS.nNodes
v2Indices=np.arange(MS.nNodes).astype(int) + 2*MS.nNodes
enIndices=np.arange(MS.nNodes).astype(int) + 3*MS.nNodes

outerIndices = []
for name in ['Inflow', 'Outflow', 'Top', 'Bottom']:
	gInd = MS.groupNames.index(name)
	outerIndices = list( set(outerIndices) | set(MS.groupMembers[gInd]) ) # logical union
outerIndices.sort()
outerIndices = np.array(outerIndices)

bodyInd = MS.groupNames.index('Body')
bodyIndices = np.array(MS.groupMembers[bodyInd]) # Don't sort these




# get Finite Element operators
GradGradBuilder = MS.makeStiffnessMatrix()
Mixed1Builder, Mixed2Builder = MS.makeMixedMatrices()

MM = GradGradBuilder.getFullSparse().tolil()
M1 = Mixed1Builder.getFullSparse().tolil()
M2 = Mixed2Builder.getFullSparse().tolil()

for i in range(NUM_VARS):
	MM[outerIndices+i*MS.nNodes,:] = 0
	M1[outerIndices+i*MS.nNodes,:] = 0
	M2[outerIndices+i*MS.nNodes,:] = 0


MM = MM.tocsr()
M1 = M1.tocsr()
M2 = M2.tocsr()



# Set initial guess
U = np.ones(NUM_VARS*MS.nNodes)
U[v1Indices] = Vinf[0]
U[v2Indices] = Vinf[1]
U[enIndices] = en_inf


# Printing inital Solution
MS.printResults(U[rhoIndices], U[v1Indices], U[v2Indices], U[enIndices], 0)



for it in range(1,numIts+1):

	# Compute flux functions at each node
	F1, F2 = FF.F12(U[rhoIndices], U[v1Indices], U[v2Indices], U[enIndices])
	# Compute Jacobian of each flux function
	DF1DU, DF2DU = FF.df12dU(U[rhoIndices], U[v1Indices], U[v2Indices], U[enIndices])

	# Compute Residual
	F = (M1 @ F1) + (M2 @ F2) - visc*( MM @ U )
	if np.linalg.norm(F) > 1e9 :
		print("Blew Up!")
		break
	
	# Compute Jacobian of Residual
	DFDU = (M1 @ DF1DU) + (M2 @ DF2DU) - visc*MM

	# Enforce Boundary Conditions
	DFDU = DFDU.tolil()
	for i in range(NUM_VARS):
		for j in range(len(outerIndices)):
			DFDU[ outerIndices[j]+i*MS.nNodes, outerIndices[j]+i*MS.nNodes ] = 1


	
	DFDU = DFDU.tocsr()

	del_U = spsolve(DFDU, -F)
	U += del_U

	# Printing Solution
	MS.printResults(U[rhoIndices], U[v1Indices], U[v2Indices], U[enIndices], it)


print(np.linalg.norm(F))
