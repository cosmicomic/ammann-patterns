import numpy as np
import matplotlib.pyplot as plt
import itertools
import math

# projection tools

def normalize(x):
    return x / np.linalg.norm(x)

def project_onto_plane(x, n):
    d = np.linalg.dot(x, n) / np.lingalg.norm(n)
    p = [d * normalize(n)[i] for i in range(len(n))]
    return [x[i] - p[i] for i in range(len(x))]

# set parameters
lambdaP = .5 * (1 + np.sqrt(5))
lambdaN = .5 * (1 - np.sqrt(5))
lsRatioP = .5 * (1 + np.sqrt(5))
lsRatioN = .5 * (1 - np.sqrt(5))
pi = np.pi

n = 5

kappa1 = 1 / (1 - lsRatioN)
kappa2 = (1 / (1 - (1 / lsRatioN)))

m1P = 1 / (1 - (1/lsRatioN)*(lsRatioP))
m2P = 1 / (1 - (lsRatioN)*(1/lsRatioP) * -1 * (1/lsRatioN))

# find Coxeter plane

C = np.array( [ [0, 0, 0, 0, 1], [1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, 1, 0] ] )

coxEigvalP = np.exp((2 * np.pi * 1j) / 5)
coxEigvalN = np.exp((-2 * np.pi * 1j) / 5)

eigvecs = np.linalg.eig(C)[1]

coxSpanP = np.array( [0.309017 - 0.951057*1j, -0.809017 - 0.587785*1j, -0.809017 + 0.587785*1j, 0.309017 + 0.951057*1j, 1] )
coxSpanN = np.array( [0.309017 + 0.951057*1j, -0.809017 + 0.587785*1j, -0.809017 - 0.587785*1j, 0.309017 - 0.951057*1j, 1] )

basis1 = normalize(coxSpanP)
basis2 = coxSpanN - (coxSpanP * (np.vdot(coxSpanP, coxSpanN.transpose()) / np.vdot(coxSpanP, coxSpanP.transpose())))
basis2 = normalize(basis2)

print(basis1, basis2)

def project_to_coxeter(basis1, basis2, vector):
    proj1 = np.vdot(vector, basis1.transpose())
    proj2 = np.vdot(vector, basis2.transpose())
    return np.array([proj1, proj2])

# project cryst roots to Coxeter plane

fRoots = list(set([root for root in itertools.permutations([-1, 1, 0, 0, 0])]))
print(fRoots)
projRoots = []

for root in fRoots:
    proj1 = np.vdot(root, basis1.transpose())
    proj2 = np.vdot(root, basis2.transpose())
    projRoots.append(np.array([proj1, proj2]))



veclist = []

def rotate(radians, vector):
    matrix = np.array([[np.cos(radians), -1 * np.sin(radians)], [np.sin(radians), np.cos(radians)]])
    return np.dot(matrix, vector)

'''x = []
y = []
for root in fRoots:
    x.append(np.dot(basis1.transpose(), root).real)
    y.append(np.dot(basis2.transpose(), root).imag)'''
#plt.scatter(x,y)
#plt.savefig('circ.png')
#plt.clf()


# make a and b vectors

rotations = []
a_0 = np.array([0,1])

for k in range(n):
    array = np.array([[np.cos((2 * pi * k) / n), -1 * np.sin((2 * pi * k) / n)], [np.sin((2 * pi * k) / n), np.cos((2 * pi * k) / n)]])
    rotations.append(array)


aP_vectors = [rotation.dot(a_0) for rotation in rotations]
bP_vectors = [-1 * (1 / lsRatioN) * a_j for a_j in aP_vectors]

all_vectors = aP_vectors + bP_vectors

x = [vector[0] for vector in all_vectors]
y = [vector[1] for vector in all_vectors]

#plt.scatter(x,y)
#plt.show()

# calculate phases

q = np.array([1,0,0,0,1])
qP = project_to_coxeter(basis1, basis2, q)
qN = project_to_coxeter(basis2, -basis1, q)

print("qP", qP)
print("qN", qN)

aN_vectors = []
bN_vectors = []

for aVector in aP_vectors:
    aN_vectors.append(np.array([aVector[1], -aVector[0]]))

for bVector in bP_vectors:
    bN_vectors.append(np.array([bVector[1], -bVector[0]]))

phases1 = []
phases2 = []

for j in range(n):
    phases1.append(np.vdot(aP_vectors[j] + lsRatioP * bP_vectors[j], qP).real)
    phases2.append(np.vdot(aN_vectors[j] + lsRatioN * bN_vectors[j], qN).real)

print(phases1, phases2)



