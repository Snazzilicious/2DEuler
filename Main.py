

from scipy.sparse import csr_matrix, lil_matrix

import MeshStructs as MS
import FluxFunctions as FF

""" Set Control Variables """
visc = 1.0
numIts = 10
numSubIts = 10

# Pick a mesh file
filename = "Mesh/2DCircle.su2"

# Load in the Mesh
MS.parseMesh(filename)


""" Set free stream conditions """
M_inf = 2.0
AoA = 10*np.pi/180.0

# don't change these
Vinf = np.array([np.cos(AoA), -np.sin(AoA)])
en_inf = ( 1.0 / ( FF.gamma*(FF.gamma-1)*M_inf*M_inf ) )

# get Finite Element operators
GradGrad = MS.makeStiffnessMatrix()
D1,D2 = MS.makeMixedMatrices()


# Shortucts to Variables
rhoIndices = np.arange(MS.nNodes).astype(int)
v1Indices=np.arange(MS.nNodes).astype(int) + MS.nNodes
v2Indices=np.arange(MS.nNodes).astype(int) + 2*MS.nNodes
enIndices=np.arange(MS.nNodes).astype(int) + 3*MS.nNodes


# Set initial guess
U = np.ones(MS.nNodes)
U[v1Indices] = Vinf[0]
U[v2Indices] = Vinf[1]
U[enIndices] = en_inf

# Printing inital Solution
MS.printResults(U[rhoIndices], U[v1Indices], U[v2Indices], U[enIndices], 0)





## Compute flux functions at each node
#F1, F2 = FF.F12(U[rhoIndices], U[v1Indices], U[v2Indices], U[enIndices])
## Compute Jacobian of each flux function
#DF1DU, DF2DU = FF.df12dU(U[rhoIndices], U[v1Indices], U[v2Indices], U[enIndices])



