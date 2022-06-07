import numpy as np


def Rmatrix(a, h, T):
    # transpose on June 7th
    result = np.array( [a] * 16 ).reshape( (4,4) )
    result[0,1] = h
    result[0,0] = -2*a -h
    result[1,1] = result[2,2] = result[3,3] = -3*a
    return result


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
