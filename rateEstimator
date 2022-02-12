import sys
import math
import numpy as np
from Bio import AlignIO
import random
import statistics

# 子孫と祖先の配列データからpairをカウント
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


# pairデータから祖先塩基をカウント
def pair2nucCount(pair):
    nucleotideCount = dict()
    #alphabet = set( "".join( pair.keys() ) )
    alphabet = set( "".join( pair.keys() ) )
    for a in alphabet:
        nucleotideCount[a] = 0
    for tag in pair.keys():
        nucleotideCount[tag[0]] += pair[tag]
    return nucleotideCount

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

# matrixのmerginalを計算し、対角成分（祖先配列）に
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

# 逆行列。
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

alignment = AlignIO.read(sys.argv[1],"fasta")
a0list, h0list = list(), list()
n = len(alignment)
for k in range(1, n):
#    print( alignment[k].id, end="\t" )
    pair = pairCount(alignment[0].seq, alignment[k].seq)
#    print( "observed difference") 
    nt =  Nt(pair)
    printMatrix(nt)
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
#    bootstrap

    a0list.append( atEstimate( x[2] ) )
    h0list.append( htEstimate( x[1], x[2] ) )

    alist, hlist = list(), list()
    bootNo = 1000 
    for i in range(bootNo):
        l = len(alignment[0])
        boot = list( random.choices( range(l), k=l ) )
        a0boot, akboot = list(), list()
        for j in boot:
            a0boot.append( alignment[0][j] )
            akboot.append( alignment[k][j] )
        pairBoot = pairCount( "".join(a0boot), "".join(akboot) )
        nt =  Nt(pairBoot)
        n0 =  N0(nt) 
        n0inv = N0inv( n0 )
        pt =  Pt( nt, n0inv )
        xboot =  wxyz ( pt )
        alist.append( atEstimate( xboot[2] ) )
        hlist.append( htEstimate( xboot[1], xboot[2]) )
#    print( atEstimate( x[2] ),statistics.pstdev( alist), sep="\t", end ="\t")
#    print( htEstimate( x[1], x[2] ),statistics.pstdev( hlist), sep="\t")

print( statistics.mean(a0list), statistics.pstdev( a0list), sep="\t", end ="\t")
print( statistics.mean(h0list), statistics.pstdev( h0list), sep="\t")
print("")

