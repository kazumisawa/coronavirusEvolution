# calculate the probability from a pair of aligned sequenes
# version 0.0 by Kazuharu Misawa 2022/06/09
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
import nucleotideProbability

def nuc2num(c):
    result = -1
    nucDict = {"c":0, "t":1, "g":2, "a":3, "C":0, "T":1, "G":2, "A":3}
    try:
        result = nucDict[ c ]
    except KeyError as e:
        pass
    return result

def sequence2vector(seq):
    result = np.zeros( ( len(seq), 4 ) )
    for i in range(len(seq)):
        pos = nuc2num( seq[i] )
        if pos >= 0: 
            result[  i, pos  ] =1
    return result

# start main #

referenceList = list()
targetList = list()

for seq_record in SeqIO.parse(sys.argv[1],"fasta"):
    referenceList.append(seq_record)

for seq_record in SeqIO.parse(sys.argv[2],"fasta"):
    targetList.append(seq_record)

a0list, h0list = list(), list()


n = len(targetList)
a = 1.48e-4
h = 1.95e-3
t = 1
P = nucleotideProbability.Pmatrix(a, h, t)

for k in range(n):
    #pairwizeAlignment
    with open("temp.fas","w") as tmpfile:
        SeqIO.write( referenceList[0], tmpfile, "fasta")
        SeqIO.write( targetList[k], tmpfile, "fasta")
    mafft_cline = MafftCommandline("mafft", input = "temp.fas")
    output1, output2 = mafft_cline()
    alignment = AlignIO.read(StringIO(output1), "fasta")
    print( alignment[0].seq, sequence2vector( alignment[0].seq) )
    print( alignment[1].seq, sequence2vector( alignment[1].seq) )
    ancestral = sequence2vector( alignment[0].seq) 
    descendant = alignment[1].seq
    logP = 0
    for i in range( np.shape(ancestral)[0] ):
        site = ancestral[i]
        pos = nuc2num( descendant[i] ) 
        if pos>=0:
            Prob =  np.dot( P, site )
            P = Prob[pos] 
            if P > 0:
                logP += math.log(P)
    print(logP)

  
