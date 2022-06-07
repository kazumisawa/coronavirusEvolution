import math
import numpy as np

def Rmatrix(a, h, t):
    # transpose on June 7th
    result = np.array( [a] * 16 ).reshape( (4,4) )
    result[0,1] = h
    result[0,0] = -2*a -h
    result[1,1] = result[2,2] = result[3,3] = -3*a
    return result


def Tmatrix(a, h, t):
    # transpose on June 7th
    result = np.array( [a] * 16 ).reshape( (4,4) )
    result[0,1] = h
    result[0,0] = -2*a -h
    result[1,1] = result[2,2] = result[3,3] = -3*a
    return result


def diffTa(a, h, t):
    result = np.zeros( (4,4) )
    b = 3*a + h
    result[1,1] = -3 * t * math.exp( -b* t )
    result[2,2] = result[3,3] = -4 * t * math.exp( -4+t )
    return result

def diffQa(a, h, t):
    result = np.array( [
                [4,     0,      0,      0],
                [2,     0,      0,      0],
                [3,     0,      0,      0],
                [3,     0,      0,      0]
        ])
    return result

def diffDa(a, h, t):
    b = 3*a + h
    return -3/(b*b)


def diffQprimea(a, h, t):
    result = np.array( [
                [0,     0,      0,      0],
                [8,     -4,     -4,     -4],
                [6,     6,      -6,     -6],
                [-3,    -3,     9,      -3]
        ])
    return result

def diffTh(a, h, t):
    result = np.zeros( (4,4) )
    b = 3*a + h
    result[1,1] = - t * math.exp( -b* t )
    return result

def diffQh(a, h, t):
    result = np.array( [
            [0, 0,  0,  0],
            [2, 0,  0,  0],
            [1, 0,  0,  0],
            [1, 0,  0,  0]
        ])
    return result

def diffDh(a, h, t):
    b = 3*a + h
    return -1/(b*b)


def diffQprimeh(a, h, t):
    result = np.array( [
                [0, 0,  0,  0],
                [4, 0,  0,  0],
                [2, 2,  -2, -2],
                [-1,    -1, 3,  -1]
        ])
    return result

a = 1.0
h = 3.0
t = 1

print( diffTa(a,h,t), diffQa(a,h,t), diffDa(a,h,t), diffQprimea(a,h,t),sep="\n" )
print( diffTh(a,h,t), diffQh(a,h,t), diffDh(a,h,t), diffQprimeh(a,h,t),sep="\n" )



nucleotide = list("CTGA")
nucDict = dict()
for i in range( len(nucleotide) ):
    nucDict[ nucleotide[i] ] = i


seq ="AATTGGCCCGG"
seqLen = len(seq)

ancestral = np.zeros( ( seqLen, len(nucleotide) ) )
descendant = np.zeros( (seqLen, len(nucleotide) ) )

for i in range(seqLen):
    c = seq[i]
    ancestral[ i, nucDict[c] ] =1

#print(ancestral)

P =  Prob(0.01,0.03,1)

print(ancestral)


for i in range(seqLen):
    site = ancestral[i]
    print( np.dot( P, site ) )
