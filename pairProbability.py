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

def sequence2vector(seq):
    nucDict = {"c":0, "t":1, "g":2, "a":3}
    result = np.zeros( ( len(seq), len(nucDict) ) )
    for i in range(len(seq)):
        try:
            result[  i, nucDict[ seq[i] ]  ] =1
        except KeyError as e:
            pass
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


a = 1.0
h = 3.0
t = 1


R = nucleotideProbability.Rmatrix(a, h, t)
T = nucleotideProbability.Tmatrix(a, h, t)
Q = nucleotideProbability.Qmatrix(a, h, t)
Qprime = nucleotideProbability.Qprimematrix(a, h, t)

diffTa = nucleotideProbability.diffTa(a, h, t)
diffQa = nucleotideProbability.diffQa(a, h, t)
diffQprimea = nucleotideProbability.diffQprimea(a, h, t)

diffTh = nucleotideProbability.diffTh(a, h, t)
diffQh = nucleotideProbability.diffQh(a, h, t)
diffQprimeh = nucleotideProbability.diffQprimeh(a, h, t)

D = nucleotideProbability.D(a, h, t)
diffDa = nucleotideProbability.diffDa(a, h, t)
diffDh = nucleotideProbability.diffDh(a, h, t)

    
print( diffTa, diffQa, diffDa, diffQprimea, sep="\n" )
print( diffTh, diffQh, diffDh, diffQprimeh, sep="\n" )



#for i in range(seqLen):
#    site = ancestral[i]
#    print( np.dot( P, site ) )

  
