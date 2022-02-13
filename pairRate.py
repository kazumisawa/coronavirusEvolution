import sys
import math
import numpy as np
import random
import re
import statistics
from io import StringIO
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align.Applications import MafftCommandline

# Count nucleotide pairs between descendant and ancestor from sequences
def pairCount(ancestor, descendant):
    alphabet = set(ancestor) | set(descendant)
    pairC = dict()
    for a in alphabet:
        for b in alphabet:
            pairC[a+b] = 0

    for i in range(len(ancestor)):
        tag = ancestor[i] + descendant[i]
        pairC[tag] = pairC[tag] + 1
    return pairC

# Count ancestral nucleotideds from pair counts
def pair2nucCount(pair):
    nucleotideCount = dict()
    #alphabet = set( "".join( pair.keys() ) )
    alphabet = set( "".join( pair.keys() ) )
    for a in alphabet:
        nucleotideCount[a] = 0
    for tag in pair.keys():
        nucleotideCount[tag[0]] += pair[tag]
    return nucleotideCount

# convert pair counts to a 4x4 matrix
def Nt(pairC):
    result = np.zeros( (4,4) )
    alphabet = "ctga"
    for i in range(len(alphabet)):
        # derived 
        a = alphabet[i]
        for j in range(len(alphabet)):
            #ancetral 
            b = alphabet[j]
            result[i][j] = pairC[b+a]
            #print(a, b, pairC[b+a])
    return result

# Calculate the merginal frequencies of the pair matrix 
# and put them to diagonal component of the ancestor matrix N0.
def N0(Nt):
    s = Nt.shape
    vector = np.zeros( s[0] )
    #↓
    for i in range(s[0]):
        #→
        for j in range(s[1]):
            vector[j] += Nt[i,j]
    result = np.zeros( s )
    for i in range( s[0] ):
        result[i,i] = vector[i]
    return result

# the inverse matrix of N0
def N0inv(N0):
    result = np.zeros_like(N0)
    for i in range( result.shape[0] ):
        result[i][i] = 1/(N0[i][i])
    return result

# P(t)を推定。祖先配列の塩基数で割る
def Pt(Nt, N0inv):
    return  Nt.dot(N0inv) 

#matrixの表示
def printMatrix(M):
    s = M.shape
    for i in range(s[0]):
        print (i, end="\t")
        for j in range(s[1]):
            print (M[i,j] , end="\t")
        print("")

# wxyxを一度に表示
def wxyz(Pt):
    S = np.matrix( [ [0,2,1,1],[1,0,1,1],[1,-1,0,1],[1,-1,1,0] ] );
    result = Pt.dot(S)
#    print("product")
#    printMatrix(result)
    return ( result[0,0]/3, result[1,1]/2, (result[2,2] + result[3,3] ) /6 )

#atの推定
def atEstimate(z):
    return -math.log( 1.0 - 4.0*z)/4

#b=(3a+h)として、btの推定
def btEstimate(xy, z):
    exp4at = 1.0 - 4.0*z 
    result = -math.log( exp4at - xy ) 
    return result

#htの推定
def htEstimate(xy, z):
    at = atEstimate(z)
    ht = btEstimate(xy, z) - 3*at
    return ht


#############start##################

referenceList = list()
targetList = list()

for seq_record in SeqIO.parse(sys.argv[1],"fasta"):
    referenceList.append(seq_record)

for seq_record in SeqIO.parse(sys.argv[2],"fasta"):
    targetList.append(seq_record)

a0list, h0list = list(), list()


n = len(targetList)
org = {"C":0, "T":1, "G":2, "A":3}
print("ID", "a", "h", sep="\t")

for k in range(n):
    #pairwizeAlignment
    with open("temp.fas","w") as tmpfile:
        SeqIO.write( referenceList[0], tmpfile, "fasta")
        SeqIO.write( targetList[k], tmpfile, "fasta")
    mafft_cline = MafftCommandline("mafft", input = "temp.fas")
    output1, output2 = mafft_cline()
    alignment = AlignIO.read(StringIO(output1), "fasta")
#    print( alignment[k].id, end="\t" )
    pair = pairCount(alignment[0].seq, alignment[1].seq)
#    print( "observed difference") 
    nt =  Nt(pair)
#    printMatrix(nt)
#    print( "ansestor again" )
    n0 =  N0(nt) 
#    printMatrix( n0 )
#    print( "n0inv" )
    n0inv = N0inv( n0 )
#    printMatrix(n0inv)
#    print( "pt" )
    pt =  Pt( nt, n0inv )
#    printMatrix(pt)
    x =  wxyz ( pt )
#    print("estimate", btEstimate( x[1], x[2] ))
    print( alignment[1].id, end="\t" )
    print( atEstimate( x[2] ), end ="\t")
    print( htEstimate( x[1], x[2] ) )

