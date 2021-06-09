# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 10:07:01 2020
@author: Sasha

READS:
    HMMER results
    AA sequences of the HMMERed ORFs
    contig expression levels
DOES:
    compares isoforms of the same contig and keeps ORFs that are more than 2 aa divergent from the closest match. When an ORf is less divergent, the more highly expressed one is kept
OUTPUTS:
    a list of retained ORF identifiers
USAGE: runs directly from the Spyder environment
"""
import sys, os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio.Alphabet import generic_dna
from Bio.Alphabet import IUPAC

taxon = 'Pygmaeoconus'
Thedir = 'E:\\Dropbox\\CONIDAE_venom_evolution\\New_Clusters'
minaaHD = 2 #This parameter defines minimal divergence between the nucleotide sequences, so that they areconsidered as coding different toxins
### functions

def HD(seq1, seq2):#seqs are of equal length
    HD = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            HD+=1
    return HD

def isSameTox(seq1, seq2, minHD):
    if len(seq1) == len(seq2) and HD(seq1, seq2) < minHD:
        return True
    elif len(seq1) == len(seq2) and HD(seq1, seq2) >= minHD:
        return False
    else:
        same = False
        if len(seq1) > len(seq2):
            ref = seq1
            query = seq2
        else:
            ref = seq2
            query = seq1
        for i in range(len(ref)-len(query)+1):
            #return ref[i:i+len(query)]
            if HD(ref[i:i+len(query)], query) < minHD:
                same = True
                break
            else:
                continue
        return same

finalORFs = open(Thedir+'\\CONIDAE_final_clusters.fa', 'w')
redlog = open(Thedir+'\\CONIDAE_ORFs_removed.txt', 'w')

### READ IN aaORFs
filelist = [thing for thing in os.listdir(Thedir) if '.fas' in thing]
for filename in filelist:
    aaORFs = {}
    clusters = {}
    counts = {}
    in1 = open(Thedir+'\\'+filename, 'r')
    for line in in1:
        line = line.strip()
        if line.startswith('>'):
            ORFID = line.split('|')[1]
            clustID = line.split('|')[0][1:]
            if ':' in line.split('|')[2]:
                count = int(line.split('|')[2].split(':')[1].split(',')[0])
            else:
                count = int(line.split('|')[2])
            if clustID in clusters.keys():
                clusters[clustID].append(ORFID)
                counts[clustID].append(count)
            else:
                clusters[clustID] = [ORFID]
                counts[clustID] = [count]
            aaORFs[ORFID] = ''
        else:
            aaORFs[ORFID] = aaORFs[ORFID] + line.replace('-', '')
    in1.close()
    print len(aaORFs)
    print clusters
    #print >>finalORFs, '>###'+filename.split('.')[0]
    for thing in clusters.keys():
        for sp in ['Pfc', 'Cll', 'Clc', 'Pyg']:
            ORFlist = []
            Seqlist = []
            explist = []
            for i in range(len(clusters[thing])):
                if sp in clusters[thing][i]:
                    ORFlist.append(clusters[thing][i])
                    Seqlist.append(aaORFs[clusters[thing][i]])
                    explist.append(counts[thing][i])
            print ORFlist
            print Seqlist
            print explist
            
            if len(ORFlist) == 0:
                toprint = []
            elif len(ORFlist) == 1:
                toprint = [True]
            else:
                toprint = [True]*len(ORFlist)
                for m in range(len(Seqlist)):
                    for n in [ind for ind in range(len(Seqlist)) if ind != m]:
                        if isSameTox(Seqlist[m], Seqlist[n], minaaHD) == True and explist[m] < explist[n]:
                            toprint[m] = False
                            print filename.split('.')[0], thing, ORFlist[m], '(', explist[m], ') --->', ORFlist[n], '(', explist[n], ')'
                            break
            print toprint
            for i in range(len(toprint)):
                if toprint[i] == True:
                    print >>finalORFs, '>'+ORFlist[i]+'\n'+aaORFs[ORFlist[i]]

finalORFs.close()
redlog.close()
