# obtain trajectory of nucleotide contents 
# version 0.3 by Kazuharu Misawa 2022/02/13
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
    # a loop for rows
    for i in range(s[0]):
        # a loop for columns
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

# P(t) is estimated by dividing 
# by the number of nucleotides in the ancestral sequence.
def Pt(Nt, N0inv):
    return  Nt.dot(N0inv) 

#print a matrix
def printMatrix(M):
    s = M.shape
    for i in range(s[0]):
        print (i, end="\t")
        for j in range(s[1]):
            print (M[i,j] , end="\t")
        print("")

# calculate w, x, y, and z.  See details in manuscript
def wxyz(Pt):
    S = np.matrix( [ [0,2,1,1],[1,0,1,1],[1,-1,0,1],[1,-1,1,0] ] );
    result = Pt.dot(S)
#    print("product")
#    printMatrix(result)
    return ( result[0,0]/3, result[1,1]/2, (result[2,2] + result[3,3] ) /6 )

#estimating the rate of non-C-to-U substitions
def atEstimate(z):
    return -math.log( 1.0 - 4.0*z)/4

#estimating b=(3a+h)
def btEstimate(xy, z):
    exp4at = 1.0 - 4.0*z 
    result = -math.log( exp4at - xy ) 
    return result

#estimating the rate of C-to-U substitions
def htEstimate(xy, z):
    at = atEstimate(z)
    ht = btEstimate(xy, z) - 3*at
    return ht

def diagonal(vector):
    n = len(vector)
    result = np.zeros( [n, n] )
    for i in range(n):
        result[i][i] = vector[i]
    return result

def Q(h, a):
    return np.matrix([ [4*a,1,0,0], [2*(a+h),-1,1,0], [3*a+h,0,0,1], [3*a+h,0,-1,-1] ])

def Qinv(h, a):
    b = 3*a+h
    tmp = np.matrix([ [1,1,1,1], [4*(b-a),-4*a,-4*a,-4*a], [2*b,2*b,-2*b,-2*b], [-b,-b,3*b,-b] ] )
    return tmp/(4*b)

def Rtime(h, a, t):
    b = 3*a+h
    A = math.exp(-4*a*t)
    B = math.exp(-b*t)
    return np.matrix( [ [1,0,0,0],[0,B,0,0],[0,0,A,0],[0,0,0,A] ] )

def Ptime(h, a, t):
    Q0 = Q(h,a)
    Qinv0 = Qinv(h,a)
    R = Rtime(h,a,t)
    return np.dot( np.dot(Q0,R), Qinv0) 

def content(seq):
    alphabet = "ctga"
    count = list()
    for c in alphabet:
        count.append( seq.count(c) )
    return count


org = {"C":0, "T":1, "G":2, "A":3}

for seq_record in SeqIO.parse(sys.argv[1],"fasta"):
    seq = re.sub("A*$","",str(seq_record.seq) )
    count = dict()
    total = 0
    for c in org:
        count[c] = seq.count(c)
        total += count[c]

nRef = diagonal( [count["C"], count["T"], count["G"], count["A"]] )

t = float(sys.argv[2])
for i in range(3,len(sys.argv)):
    with open(sys.argv[i]) as f:
        for line in f:
            dollar = line.split("\t")
            ID = dollar[0]
            if ID=="ID":
                print("ID","eTime", "C", "T", "G", "A", sep="\t")
            else:
                a = float(dollar[1])
                h = float(dollar[2])
                P = Ptime(h,a,t)
                result =  merginal( np.dot(P,nRef).T ) 
                total =  result[0] + result[1] + result[2] + result[3];
                print( ID, t, result[0]/total, result[1]/total, result[2]/total, result[3]/total, sep="\t" )

